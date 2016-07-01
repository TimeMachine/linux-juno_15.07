#include "sched.h"
#include <linux/cpufreq.h>
#include <linux/slab.h>
#include <linux/sort.h>
#include <linux/hrtimer.h>
#include <linux/ktime.h>
#include <linux/kthread.h>
#include <linux/fs.h>
#include <asm/segment.h>
#include <asm/uaccess.h>
#include <linux/buffer_head.h>
////#define _debug
//#define _debug_pick
//#define _sched_debug
#define kHZ 1000
#define juno

#ifdef juno
#define _1s_jiffies 50000000
#define little_cpu 4
#define big_cpu 2
#define _MAX_CPU 6
static int index[2][4] = {{0,3,4,5}, {1,2}}; // Little core cluster and Big core cluster.
#endif

#define jiffies_to_cycle(jiffies,freq)  (((u64) jiffies * freq) / (_1s_jiffies / kHZ))
#define cycle_to_jiffies(cycle,freq)  (((u64) cycle * _1s_jiffies) / ((u64) freq * kHZ))

// lock by main_schedule
static DEFINE_SPINLOCK(mr_lock);
static spinlock_t rq_lock[_MAX_CPU];
static struct hrtimer hr_timer;
static int total_job = 0;

extern unsigned int get_stats_table(int cpu, unsigned int **freq);
//extern void change_governor_userspace(int cpu);
extern struct cpufreq_policy *cpufreq_cpu_get(unsigned int cpu);
static void main_schedule(int workload_predict);
static void put_prev_task_energy(struct rq *rq, struct task_struct *prev);

int is_big_cluster(int cpu){
	return cpu == 1 || cpu == 2;
}

static inline void _all_rq_lock(void)
{
	int i = 0;
	spin_lock(&rq_lock[i]);	
	for (i = 1; i < _MAX_CPU; i++) 
		spin_lock_nested(&rq_lock[i], i);	
}

static inline void _all_rq_unlock(void)
{
	int i;
	for (i = 0; i < _MAX_CPU; i++) 
		spin_unlock(&rq_lock[i]);	
}

static inline void _double_rq_lock(int lock1, int lock2)
{
	if (lock1 == lock2)	
		spin_lock(&rq_lock[lock1]);	
	else if (lock1 < lock2) {
		spin_lock(&rq_lock[lock1]);	
		spin_lock_nested(&rq_lock[lock2], SINGLE_DEPTH_NESTING);	
	}
	else {
		spin_lock(&rq_lock[lock2]);	
		spin_lock_nested(&rq_lock[lock1], SINGLE_DEPTH_NESTING);
	}
}

static inline void _double_rq_unlock(int lock1, int lock2)
{
	spin_unlock(&rq_lock[lock1]);	
	if (lock1 != lock2)	
		spin_unlock(&rq_lock[lock2]);	
}

static unsigned long long get_cpu_cycle(void)
{
	return get_cycles();
	//return cpu_cycle();
}

static void get_cpu_frequency(int cpu)
{
	struct energy_rq *e_rq = &cpu_rq(cpu)->energy;
	struct cpufreq_policy *policy = cpufreq_cpu_get(cpu);		
	/*
	struct cpufreq_stats *stats = policy->stats;
	BUG_ON(policy == NULL);

	e_rq->freq = stats->freq_table;
	e_rq->state_number = stats->state_num;		
	*/
	e_rq->state_number = get_stats_table(cpu, &e_rq->freq);
	e_rq->freq_now = policy->cur;
	cpufreq_cpu_put(policy);
	/*int i = 0;	
	for (i = 0; i < e_rq->state_number; i++) {
		printk("cpu:%d  freq[%d] %u\n", cpu, i, e_rq->freq[i]);
	}*/
}

static void reschedule(int cpu)
{
	struct rq *i_rq = cpu_rq(cpu);
	//printk("resched %d, pid:%d\n",cpu,i_rq->curr->pid);
	if (cpu == smp_processor_id())
		set_tsk_need_resched(i_rq->curr);
	else if(cpu != smp_processor_id() && i_rq->curr == i_rq->idle) 
		// test....
		wake_up_nohz_cpu(cpu);
		//wake_up_idle_cpu(cpu);
	else
		resched_cpu(cpu);
}

static enum hrtimer_restart sched_period_timer(struct hrtimer *timer)
{
	ktime_t ktime;
	//get_cycles()
	//printk("[period] cycle:%llu total_job:%d\n",get_cpu_cycle(), total_job);

	ktime = ktime_set(0, NSEC_PER_SEC);

	spin_lock(&mr_lock);
	main_schedule(true);
	spin_unlock(&mr_lock);
	
	hrtimer_forward_now(timer, ktime);
	
	return HRTIMER_RESTART;
}

static void update_credit(struct task_struct *curr)
{
	if (curr->ee.execute_start != 0) {
		int cpu = task_cpu(curr);
		u64 executionTime = get_cpu_cycle() - curr->ee.execute_start;
		curr->ee.total_execution += executionTime;
		if (curr->ee.credit[cpu] > executionTime)
 			curr->ee.credit[cpu] -= executionTime; 
		else
			curr->ee.credit[cpu] = 0;
	}		
	curr->ee.execute_start = get_cpu_cycle();
}

static void update_curr_task(struct rq *rq)
{
	struct task_struct *curr = rq->curr;
	u64 delta_exec;
	delta_exec = rq->clock_task - curr->se.exec_start;
	if (unlikely((s64)delta_exec < 0))
		delta_exec = 0;

	schedstat_set(curr->se.statistics.exec_max,
			max(curr->se.statistics.exec_max, delta_exec));

	curr->se.sum_exec_runtime += delta_exec;
	account_group_exec_runtime(curr, delta_exec);

	curr->se.exec_start = rq->clock_task;
	cpuacct_charge(curr, delta_exec);
}

static void update_curr_energy(struct rq *rq)
{
	
	//update the energy info.
	update_credit(rq->curr);
	rq->energy.time_sharing = rq->clock_task;
	//update the curr info.		
	update_curr_task(rq);
}

static inline int on_energy_rq(struct sched_energy_entity *ee)
{
    return !list_empty(&ee->list_item);
}

static void move_task_to_rq(struct rq *rq, struct sched_energy_entity *task ,int set)
{
	if (task_contributes_to_load(task->instance))
		task->rq_e->rq->nr_uninterruptible++;
	list_del(&task->list_item);
	task->rq_e->energy_nr_running--;
	dec_nr_running(task->rq_e->rq);
	
	//printk("before move pid:%d energy_nr:%lu nr:%u task_contributes_to_load:%d cpu:%d\n",task->instance->pid, task->rq_e->energy_nr_running, task->rq_e->rq->nr_running, task_contributes_to_load(task->instance),task->rq_e->rq->cpu);
	if (set)
		set_task_cpu(task->instance, rq->cpu);
	
	if (task_contributes_to_load(task->instance))
		rq->nr_uninterruptible--;
	list_add_tail(&task->list_item,&rq->energy.queue);
	task->rq_e = &rq->energy;
	task->rq_e->energy_nr_running++;
	inc_nr_running(rq);
	
	//printk("after move pid:%d energy_nr:%lu nr:%u task_contributes_to_load:%d cpu:%d\n",task->instance->pid, task->rq_e->energy_nr_running, task->rq_e->rq->nr_running, task_contributes_to_load(task->instance),task->rq_e->rq->cpu);
}

static int _migration(void *data)
{
	struct migration_arg *arg = data;
	set_task_cpu(arg->task, arg->dest_cpu);
	//move_task_to_rq(cpu_rq(arg->dest_cpu), &arg->task->ee, 1);
	return 0;
}

static void rq_selection(int *try_cpu, int this_cpu)
{
	int cluster = 0;
	int i;
	int core_count;
	for (i = 0; i < big_cpu; i++)
		if (this_cpu == index[1][i])	
			cluster = 1;
	core_count = cluster == 0 ? little_cpu : big_cpu;
	/*
	 * Priority : split task -> non-split task -> other CPU run queue
	 * try_cpu  :		0				1			    2 , 3
	 */
	try_cpu[0] = this_cpu;
	try_cpu[1] = this_cpu;
	// try to steal other run queue
	try_cpu[2] = index[cluster][(this_cpu - 1 + core_count) % core_count];
	try_cpu[3] = index[cluster][(this_cpu + 1) % core_count];
}

#ifdef CONFIG_SMP
static int
select_task_rq_energy(struct task_struct *p, int sd_flag, int flags)
{
	return task_cpu(p); 
}

static void post_schedule_energy(struct rq *rq)
{
    rq->curr->ee.select = 0;
}

static void pre_schedule_energy(struct rq *rq, struct task_struct *prev)
{
#ifdef _debug
	//printk("cpu:%d, %s ,pid:%d ,need_move:%d\n",smp_processor_id(),__PRETTY_FUNCTION__,prev->pid,prev->ee.need_move);
#endif
	prev->ee.select = 0;
			
	if (prev->ee.need_move == -1 && 
		prev->on_rq &&
		rq->cpu == task_cpu(prev) && 
		prev->ee.split && 
		prev->ee.credit[rq->cpu] == 0) {
		int move = -1;
		int near[4];
		int i;
		_all_rq_lock();
		rq_selection(near, rq->cpu);
		for (i = 2; i <= 3; i++) {
			if (prev->ee.credit[near[i]])
				move = near[i];
		}
		if (move >= 0) {
#ifdef _debug_pick
			printk("[need move] :%d\n",move);
#endif
			move_task_to_rq(cpu_rq(move) ,&prev->ee ,0);
			reschedule(move);
		}
		_all_rq_unlock();
	}
}

static void task_waking_energy(struct task_struct *p)
{
	struct rq *rq = p->ee.rq_e->rq;
	p->ee.credit[rq->cpu] = cycle_to_jiffies(p->ee.workload_guarantee * kHZ, rq->energy.freq_now);
	//p->ee.total_execution = 0;
}

#endif /* CONFIG_SMP */

static void
check_preempt_curr_energy(struct rq *rq, struct task_struct *p, int flags)
{
	/* we're never preempted */
}

static struct task_struct *
pick_next_task_energy(struct rq *rq)
{
	struct task_struct *next;
	struct sched_energy_entity *next_ee;
	struct energy_rq *e_rq;
	struct list_head *pos;
	int find = 0, i = 0, j = 0;
	int try_cpu[4] = {0};
	int next_ee_cpu = 0;
	struct sched_energy_entity *retry = NULL;
	rq_selection(try_cpu, rq->cpu);	

	for (i = 0; i < 4; i++) {
		_all_rq_lock();
		e_rq = &cpu_rq(try_cpu[i])->energy;
		//if(cpu_rq(try_cpu[i])->online != 1 || e_rq->energy_nr_running == 0){
		if(e_rq->energy_nr_running == 0){
			//if (rq->curr->policy == SCHED_ENERGY)
			//	printk("[debug pick energy_nr_running=0] cpu:%d try:%d \n",e_rq->rq->cpu,try_cpu[i]);
			_all_rq_unlock();
			continue;
		}
		pos = e_rq->queue.next;
		for (j = 0; pos != &e_rq->queue; pos = pos->next) {
			next_ee = list_entry(pos, struct sched_energy_entity, list_item);
			next_ee_cpu = task_cpu(next_ee->instance);
			//if (i == 0||i==1)
			//printk("[debug pick] cpu:%d pid:%d i:%d nr:%lu %d %d %d %d %d\n",rq->cpu,next_ee->instance->pid,i,e_rq->energy_nr_running,!next_ee->credit[rq->cpu],next_ee->instance->state != TASK_RUNNING,next_ee->select == 1,next_ee->need_move >= 0 && next_ee->need_move != rq->cpu,next_ee_cpu != try_cpu[i] &&(cpu_rq(next_ee_cpu)->curr->pid == next_ee->instance->pid));
			if (!next_ee->credit[rq->cpu])
				continue;

			if (next_ee->instance->state != TASK_RUNNING)
				continue;

			if (next_ee->select == 1)
				continue;

			if (next_ee->need_move >= 0 && next_ee->need_move != rq->cpu)
				continue;

			if (next_ee_cpu != try_cpu[i] && 
					(cpu_rq(next_ee_cpu)->curr->pid == next_ee->instance->pid)) { 
				if (i == 1) {
					//printk("[pick fail] pid:%d cpu:%d next_ee_cpu:%d try_cpu:%d curr_pid:%d locked:%d\n", next_ee->instance->pid,rq->cpu,next_ee_cpu,try_cpu[i],cpu_rq(next_ee_cpu)->curr->pid,raw_spin_is_locked(&cpu_rq(next_ee_cpu)->lock));
					// if we not found any task, we will yield it to pick.
					retry = next_ee;
				}
				continue;
			}

			if (i <= 1) {
				if ((i == 0 && next_ee->split) || i == 1) {
					find = 1;
					//if (next_ee_cpu != rq->cpu) {
						set_task_cpu(next_ee->instance, rq->cpu);
					//}
					break;
				}
			}
			//else if (!raw_spin_is_locked(&cpu_rq(try_cpu[i])->lock) && 
			else if (next_ee->instance->pid != cpu_rq(try_cpu[i])->curr->pid) {
				//printk("steal | i:%d, next-pid:%d, curr-pid:%d, j:%d\n",try_cpu[i],next_ee->instance->pid, cpu_rq(try_cpu[i])->curr->pid,j);
				//steal other cpu run queue.
				move_task_to_rq(rq ,next_ee ,1);
				find = 1;
				break;
			}
			//printk("[picking] cpu:%d try:%d pid:%d loop:%d need_move:%d state:%ld\n",rq->cpu,try_cpu[i],next_ee->instance->pid,j++,next_ee->need_move ,next_ee->instance->state);
		}

		if(find == 1) {
			// To be executed task has put first entry.
			next_ee->select = 1;
			next_ee->need_move = -1;
			_all_rq_unlock();
			break;
		}
		_all_rq_unlock();
	}
	if (find == 0) {
#ifdef _debug_pick
		if (rq->curr->policy == SCHED_ENERGY)
			printk("[pick NULL] cpu:%d \n",smp_processor_id());
#endif
		/*if (retry != NULL) {
			// yield the task and pick.
			migrate_task(retry->instance, rq->cpu);
			//	reschedule(task_cpu(retry->instance));
			//set_task_cpu(retry->instance, rq->cpu);
			_all_rq_lock();
			next_ee = retry;
			next_ee->select = 1;
			next_ee->need_move = -1;
			_all_rq_unlock();
		}
		else*/
			return NULL;
	}
	next = container_of(next_ee, struct task_struct, ee);
	next->se.exec_start = rq->clock_task;
	next->ee.execute_start = get_cpu_cycle();
	rq->post_schedule = 1;
#ifdef _debug_pick
	printk("[%s] cpu:%d pid:%d credit:%llu\n",__PRETTY_FUNCTION__,rq->cpu,next->pid, next->ee.credit[rq->cpu]);
#endif
	return next;
}

static void
enqueue_task_energy(struct rq *rq, struct task_struct *p, int flags)
{
	int i = 0;
	int cpu = rq->cpu;
	struct rq *true_rq;
#ifdef _debug
	printk("cpu:%d, %s ,pid:%d nr:%lu begin\n",smp_processor_id(),__PRETTY_FUNCTION__,p->pid,rq->energy.energy_nr_running);
#endif
	//printk("[debug nr] nr:%d fair_nr:%d\n",rq->nr_running, rq->cfs.h_nr_running);
	// for the first time (init cpu freq)
	if ( unlikely(rq->energy.freq == NULL) ) {
		for (i = 0; i < _MAX_CPU; i++) {  
			get_cpu_frequency(i);
		}
	}
	_all_rq_lock();
	if (p->ee.need_move != -1)
		cpu = p->ee.need_move;
	true_rq = cpu_rq(cpu);

	list_add_tail(&(p->ee.list_item),&(true_rq->energy.queue));
	p->ee.rq_e = &true_rq->energy;
	true_rq->energy.energy_nr_running++;
	inc_nr_running(true_rq);
	//printk("*enqueue* pid:%d cpu:%d nr:%lu flag:%d\n",p->pid,cpu,true_rq->energy.energy_nr_running, flags);

	if((flags == 0 || flags == 2) && p->ee.first == 1) {
	//update the new task info.
		if (p->ee.workload_guarantee)
			p->ee.workload = p->ee.workload_guarantee * kHZ;
		else
			p->ee.workload = rq->energy.freq[0] * kHZ;
		p->ee.credit[rq->cpu] = cycle_to_jiffies(p->ee.workload, rq->energy.freq_now);
		//printk("[enqueue] new job pid:%d, total_job:%d\n",p->pid, total_job);
		total_job++;
		if (total_job == 1) {
			ktime_t ktime;
			ktime = ktime_set(0, NSEC_PER_SEC);
			start_bandwidth_timer(&hr_timer, ktime);
		}
		p->ee.first = 0;
	}
	_all_rq_unlock();
#ifdef _debug
	printk("[debug en] pid:%d prev->state:%ld flag:%d\n", p->pid, p->state, flags);
#endif
}

static void
dequeue_task_energy(struct rq *rq, struct task_struct *p, int flags)
{
	int cpu = rq->cpu;
	struct rq *true_rq;
	_all_rq_lock();
#ifdef _debug
	printk("cpu:%d, %s ,pid:%d nr:%lu need_move:%d begin\n",smp_processor_id(),__PRETTY_FUNCTION__,p->pid,rq->energy.energy_nr_running,p->ee.need_move);
#endif
	//printk("[debug nr] nr:%d fair_nr:%d\n",rq->nr_running, rq->cfs.h_nr_running);
	if (p->ee.need_move != -1)
		cpu = p->ee.need_move;
	true_rq = cpu_rq(cpu);

	update_curr_energy(true_rq);

	list_del(&(p->ee.list_item));
	true_rq->energy.energy_nr_running--;
	dec_nr_running(true_rq);
	p->ee.select = 0;
	_all_rq_unlock();

	if(flags) {
		if (p->state == TASK_DEAD) {
			total_job--;
			if (total_job == 0) {
				hrtimer_cancel(&hr_timer);
			}
		}
	}
#ifdef _debug
	printk("[debug de] pid:%d prev->state:%ld flag:%d\n", p->pid, p->state, flags);
#endif
}

static void requeue_task_energy(struct rq *rq, struct task_struct *p)
{
	BUG_ON(rq->cpu != p->ee.rq_e->rq->cpu);
	if (on_energy_rq(&p->ee)) {
		list_move_tail(&p->ee.list_item, &rq->energy.queue);
	}
}

static void yield_task_energy(struct rq *rq)
{
	update_curr_energy(rq);
}

static void put_prev_task_energy(struct rq *rq, struct task_struct *prev)
{
#ifdef _debug
	//printk("cpu:%d, %s\n",smp_processor_id() ,__PRETTY_FUNCTION__);
#endif
	update_curr_energy(rq);
}

static unsigned int next_frequency_level(unsigned int workload, int cpu)
{
	struct energy_rq *i_rq = &cpu_rq(cpu)->energy;
	int i;
	for (i = 0; i < i_rq->state_number; i++) {
		if (i_rq->freq[i] * kHZ > workload)
			return i_rq->freq[i] * kHZ;
	}
	if (cpu == 1 || cpu == 2)
		return i_rq->freq[i-1] * kHZ;
	if (workload > 850000 * kHZ && workload <= 950000 * kHZ)
		return 950000 * kHZ;
	return 1100000 * kHZ;
}

static void workload_prediction(void)
{
	struct rq *i_rq;
	int i = 0;
	struct list_head *head;
	struct sched_energy_entity *data;
	struct list_head *pos;
	int next_level; 

	_all_rq_lock();
	for (i = 0 ;i < _MAX_CPU; i++) {
		i_rq = cpu_rq(i);
		if (i_rq->energy.energy_nr_running != 0) {
			head = &i_rq->energy.queue;
			for(pos = head->next; pos != head; pos = pos->next) {
				data = list_entry(pos ,struct sched_energy_entity, list_item);
				next_level = next_frequency_level(data->workload, i);
				// predict workload from the statistics
				//printk("pid:%d, cpu:%d, exeute_start:%llu, total_execution:%llu, workload:%llu, workload(jiffies):%llu, freq_now:%u\n",data->instance->pid , i, data->execute_start ,data->total_execution, data->workload, cycle_to_jiffies(data->workload, i_rq->energy.freq_now), i_rq->energy.freq_now);
				if (data->workload > i_rq->energy.freq_now * kHZ) {
					// special case (cpufreq is not setting).
					data->workload = jiffies_to_cycle(data->total_execution, i_rq->energy.freq_now);
					data->over_predict = 0;
				}
				else if (data->total_execution == 0) {
					// for the newly job
					if (data->workload_guarantee)
						data->workload = data->workload_guarantee * kHZ;
					else
						data->workload = i_rq->energy.freq[0] * kHZ;
					data->over_predict = 0;
				}
				else if (data->total_execution >= cycle_to_jiffies(data->workload, i_rq->energy.freq_now) * 9 / 10) {
				//else if (jiffies_to_cycle(data->total_execution, i_rq->energy.freq_now) >= data->workload) {
					if (data->workload_guarantee)
						data->workload = data->workload_guarantee * kHZ;
					else {
						if (data->over_predict) 
							data->workload = next_level;
						//data->workload = (i_rq->energy.freq[i_rq->energy.state_number - 1] + data->workload / kHZ) / 2 * kHZ;
						else {
							int increase = 0;
							if (i == 1||i == 2)
								increase = 175000 * kHZ;
							else
								increase = 125000 * kHZ;
							if (data->workload + increase >= next_level)
								data->workload = next_level;
							else 
								data->workload += increase;  
							//if (data->workload > i_rq->energy.freq[i_rq->energy.state_number - 1])
							//	data->workload = i_rq->energy.freq[i_rq->energy.state_number - 1] *kHZ;
						}
					}
					data->over_predict++;
				}
				else {
					//if (is_big_cluster(i))
					//	data->workload = (jiffies_to_cycle(data->total_execution, i_rq->energy.freq_now) + (data->workload)*3) / 4;
					//else
						data->workload = ((jiffies_to_cycle(data->total_execution, i_rq->energy.freq_now) + data->workload)) / 2;
					data->over_predict = 0;
				}
				// reset the statistics.
				data->total_execution = 0;
			}
		}
	}

	_all_rq_unlock();
}

static int compare_alpha(const void *a, const void *b)
{
	struct sched_energy_entity * const *ippa = a;	
	struct sched_energy_entity * const *ippb = b;	
	const struct sched_energy_entity *ipa = *ippa;	
	const struct sched_energy_entity *ipb = *ippb;	
	if ((ipa->dummy_workload > (u64)850000 * 10000 / 17) && is_big_cluster(task_cpu(ipa->instance))
			&& (ipb->dummy_workload <= 85000 * kHZ) && !is_big_cluster(task_cpu(ipb->instance)) )
		return -1;
	if ((ipa->dummy_workload <= 850000 * kHZ) && !is_big_cluster(task_cpu(ipb->instance))
			&& (ipb->dummy_workload > (u64)85000 * 10000 / 17) && is_big_cluster(task_cpu(ipa->instance)))
		return 1;
	/*if (ipa->alpha > ipb->alpha)
		return -1;
	if (ipa->alpha < ipb->alpha)
		return 1;*/
	if (ipa->dummy_workload > ipb->dummy_workload)
		return -1;
	if (ipa->dummy_workload < ipb->dummy_workload)
		return 1;
	return 0;
}

static int compare_workload(const void *a, const void *b)
{
	struct sched_energy_entity * const *ippa = a;	
	struct sched_energy_entity * const *ippb = b;	
	const struct sched_energy_entity *ipa = *ippa;	
	const struct sched_energy_entity *ipb = *ippb;	
	if (ipa->dummy_workload > ipb->dummy_workload)
		return -1;
	if (ipa->dummy_workload < ipb->dummy_workload)
		return 1;
	return 0;
}
#define Max_jobs 1000
static struct sched_energy_entity *data[Max_jobs]; 
static unsigned int o_freq[_MAX_CPU];
static void scheduling_cluster(int *cpu_mask, int core_count, int job_count, int ptr,
								 u64 total_workload, unsigned int *o_freq)
{
	#define soft_float 1000000
	int i = 0, j = 0, k = 0;
	int ptr_max = ptr;
	unsigned int pre_load = 0;
	struct rq *i_rq;
	int base = ptr;

	pre_load = 0;
	for (; k < core_count; k++) {
		int f_total = 0;
		int a_jp = 0;
		i = cpu_mask[k];
		//printk("[algo debug] i:%d ptr:%d total_workload:%llu\n",i ,ptr,total_workload);
		i_rq = cpu_rq(i);
		if (total_workload == 0) {
			o_freq[i] = 0;
			continue;
		}
		for (j = i_rq->energy.state_number - 1; j >= 0; j--) {
			if (((u64) i_rq->energy.freq[j] * kHZ < data[ptr_max]->dummy_workload) ||
					((u64) i_rq->energy.freq[j] * kHZ * (core_count - k) < total_workload) ||
					((1 * soft_float - pre_load) < (data[ptr]->dummy_workload) / (i_rq->energy.freq[j] * kHZ / soft_float) ))
				break;
		}
		//printk("__debug cpu_rq(i)->energy.freq[j]:%d,data[ptr_max]->workload:%llu,total_workload:%llu,pre_load:%u\n",cpu_rq(i)->energy.freq[j],data[ptr_max]->dummy_workload,total_workload,pre_load);
		o_freq[i] = j == 4 ? i_rq->energy.freq[4] : i_rq->energy.freq[j + 1];  	
		f_total = o_freq[i] * kHZ;
		a_jp = 0;
		while ( ptr < base + job_count ) {
			int state = 0;
			if (f_total >= data[ptr]->dummy_workload)
				state = 1;
			a_jp = state ? data[ptr]->dummy_workload : f_total;
			total_workload -= a_jp;
			f_total -= a_jp;
			data[ptr]->credit[i] = cycle_to_jiffies(a_jp, o_freq[i]);
			if (state) {
				if (data[ptr]->need_move != -1 && data[ptr]->need_move != i) {
					data[ptr]->need_move = i;
					move_task_to_rq(i_rq ,data[ptr] ,0);
				}
				else {
					if (data[ptr]->rq_e->rq->cpu != i) {
					//if (data[ptr]->rq_e->rq->cpu != i && data[ptr]->instance->on_rq) {
						// if task is scheduled to other cpu, task have to move to other queue.
						//printk("[algo] pid:%d ,from:%d ,to:%d\n",data[ptr]->instance->pid,data[ptr]->rq_e->rq->cpu,i_rq->cpu);
						if (data[ptr]->select == 0 && data[ptr]->instance->state == TASK_RUNNING && 
								data[ptr]->instance->pid != data[ptr]->rq_e->rq->curr->pid &&
								data[ptr]->instance->pid != cpu_rq(task_cpu(data[ptr]->instance))->curr->pid)
							move_task_to_rq(i_rq ,data[ptr] ,1);
						else {
							data[ptr]->need_move = i;
							move_task_to_rq(i_rq ,data[ptr] ,0);
						}
						//printk("[algo] select:%d from_pid:%d true_pid:%d\n",data[ptr]->select,data[ptr]->rq_e->rq->curr->pid,cpu_rq(task_cpu(data[ptr]->instance))->curr->pid);
					}
					//else
					//	printk("[algo] pid:%d cpu:%d i:%d on_rq:%d\n",data[ptr]->instance->pid,data[ptr]->rq_e->rq->cpu,i,data[ptr]->instance->on_rq);
				}
				ptr++;
				ptr_max = ptr;
				pre_load = 0;
			}
			else {
				data[ptr]->dummy_workload -= a_jp;
				data[ptr]->split = 1;
				pre_load = a_jp / (o_freq[i] * kHZ / soft_float);
				if (ptr + 1 < base + job_count && 
						data[ptr]->dummy_workload < data[ptr+1]->dummy_workload)
					ptr_max++;
			}
			if (!f_total)
				break;
		}
	}
}

static void algo(int workload_predict)
{
	u64 total_workload = 0;
	int job_count = 0, cluster_job = 0;
	struct rq *i_rq;
	struct list_head *head;
	struct list_head *pos;
	int i = 0 ,j = 0, k = 0;
	int core_count = 0;

	_all_rq_lock();
	for (i = 0 ;i < _MAX_CPU; i++) {
		i_rq = cpu_rq(i);
		if (i_rq->energy.energy_nr_running != 0) {
			head = &i_rq->energy.queue;
			for(pos = head->next; pos != head; pos = pos->next) {
				data[job_count] = list_entry(pos ,struct sched_energy_entity, list_item);
				if (data[job_count]->instance->state == TASK_RUNNING && 
					data[job_count]->workload > jiffies_to_cycle(data[job_count]->total_execution, i_rq->energy.freq_now)) {
					if (workload_predict)
						data[job_count]->dummy_workload = data[job_count]->workload;
					else
						data[job_count]->dummy_workload = data[job_count]->workload - jiffies_to_cycle( data[job_count]->total_execution, i_rq->energy.freq_now);
					for (k = 0 ;k < _MAX_CPU; k++) 
						data[job_count]->credit[k] = 0;
					data[job_count]->split = 0;
					job_count++;
				}
			}
		}
	}
	if (job_count == 0) {
		_all_rq_unlock();
		return;
	}
	memset(o_freq, 0, _MAX_CPU * sizeof(unsigned int));
	sort(data, job_count, sizeof(struct sched_energy_entity*), compare_alpha, NULL);
#ifdef _sched_debug
	printk("========= input===========\n");
	for(i=0;i<_MAX_CPU;i++){
		printk("CPU %d:", i);
		for(j=0;j<cpu_rq(i)->energy.state_number;j++)
			printk(" %4d", cpu_rq(i)->energy.freq[j]);
		printk("\n");
	}
	printk("job workload:");
	for(j=0;j<job_count;j++)
		printk("  (%d) %4llu %d|", data[j]->instance->pid, data[j]->dummy_workload, data[j]->alpha);
	printk("\n");
#endif

	for (j = 0, k = 0; k < 2 && j != -1; k++) {
		core_count = k == 0 ? little_cpu : big_cpu;
		j = job_count - cluster_job - 1; //the index of the remaining task.
		cluster_job = 0;
		total_workload = 0;
		for (; j >= 0; j--) {
			//over big cluster workload. Given tasks as more as possible.
			if (total_workload + data[j]->dummy_workload > 
				(u64)cpu_rq(index[k][0])->energy.freq[cpu_rq(index[k][0])->energy.state_number - 1] * core_count * kHZ && k == 0) 
				break;
			if (k == 0 && ((data[j]->dummy_workload > 850000 * kHZ && !is_big_cluster(task_cpu(data[j]->instance)))
						|| (data[j]->dummy_workload > (u64)850000 * 10000 / 17 && is_big_cluster(task_cpu(data[j]->instance)))))
				break;
			if (k == 0 && is_big_cluster(task_cpu(data[j]->instance)) 
				&& data[j]->workload > (u64)850000 * 10000 / 17) {
				data[j]->workload = 850000 * kHZ;
				data[j]->dummy_workload = data[j]->workload; 
			}
			if (k == 1 && !is_big_cluster(task_cpu(data[j]->instance))
				&& data[j]->dummy_workload > 850000 * kHZ) {
				data[j]->workload = 800000 * kHZ;
				data[j]->dummy_workload = data[j]->workload; 
			}
			total_workload += data[j]->dummy_workload;
			cluster_job++;
		}
		if (total_workload == 0)
			continue;
		//printk("[debug] k:%d j:%d cluster_job:%d total_workload:%llu\n",k,j+1,cluster_job,total_workload);
		sort(data + j + 1, cluster_job, sizeof(struct sched_energy_entity*), compare_workload, NULL);
		scheduling_cluster(index[k], core_count, cluster_job, j+1, total_workload, o_freq);	
	}
#ifdef _sched_debug
	printk("=========output===========\n");
	for(i=0;i<_MAX_CPU;i++)
		printk("CPU:%d freq:%d\n",i , o_freq[i]);
	for(i=0;i<job_count;i++){
		printk("job:%d, credit:", i);
		for(j=0;j<_MAX_CPU;j++)
			printk(" %4llu", data[i]->credit[j]);
		printk("\n");
	}
#endif

	// symmetric in the cluster.
	for(i = 0; i < 2; i++) { 
		int max_freq = -1;
		core_count = i == 0 ? little_cpu : big_cpu;
		for (j = 0; j < core_count; j++)
			if (max_freq < (int)o_freq[index[i][j]])
				max_freq = o_freq[index[i][j]];
		if (max_freq == 0)
            max_freq = cpu_rq(i)->energy.freq[0];
		if (max_freq != cpu_rq(i)->energy.freq_now)
			cpu_rq(i)->energy.set_freq = max_freq;	
	}
	
	for (i = 0 ;i < _MAX_CPU; i++) {
		i_rq = cpu_rq(i);
		if (i_rq->curr->ee.credit[i] == 0 || i_rq->curr->ee.need_move != -1){
			reschedule(i);
		}
	}
	_all_rq_unlock();
}

static void main_schedule(int workload_predict)
{
	
	struct rq *i_rq;
	int i = 0;
	// update all job data and then use them to predict.
	for (i = 0 ;i < _MAX_CPU; i++) {
		i_rq = cpu_rq(i);
		if (i_rq->energy.energy_nr_running != 0)
			update_curr_energy(i_rq);
	}	
	if (total_job) {
		if (workload_predict == true) 
			workload_prediction();
		algo(workload_predict);
	}
}

static void task_tick_energy(struct rq *rq, struct task_struct *curr, int queued)
{
#ifdef _debug
		//printk("cpu:%d, %s, pid:%d\n",smp_processor_id() ,__PRETTY_FUNCTION__, curr->pid);
#endif
	update_credit(curr);
	update_curr_task(rq);
	if (!spin_is_locked(&mr_lock) && 
		((rq->clock_task - rq->energy.time_sharing >= 30 * USEC_PER_SEC)||
		curr->ee.credit[rq->cpu] == 0)){
		// time sharing
		rq->energy.time_sharing = rq->clock_task;
		spin_lock(&rq_lock[rq->cpu]);	
		if (curr->ee.need_move == -1) 
			requeue_task_energy(rq,curr);
		spin_unlock(&rq_lock[rq->cpu]);	
		reschedule(rq->cpu);
	}
}

static void set_curr_task_energy(struct rq *rq)
{
	struct task_struct *p = rq->curr;

	p->se.exec_start = rq->clock_task;
}

static void switched_to_energy(struct rq *rq, struct task_struct *p)
{
}

static void
prio_changed_energy(struct rq *rq, struct task_struct *p, int oldprio)
{
}

static unsigned int
get_rr_interval_energy(struct rq *rq, struct task_struct *task)
{
	return 0;
}

const struct sched_class energy_sched_class = {
	.next			= &rt_sched_class,

	.enqueue_task		= enqueue_task_energy,
	.dequeue_task		= dequeue_task_energy,
	.yield_task		= yield_task_energy,

	.check_preempt_curr	= check_preempt_curr_energy,

	.pick_next_task		= pick_next_task_energy,
	.put_prev_task		= put_prev_task_energy,

#ifdef CONFIG_SMP
	.select_task_rq		= select_task_rq_energy,
	.pre_schedule       = pre_schedule_energy,
	.post_schedule      = post_schedule_energy,
	.task_waking        = task_waking_energy,
#endif

	.set_curr_task          = set_curr_task_energy,
	.task_tick		= task_tick_energy,

	.get_rr_interval	= get_rr_interval_energy,

	.prio_changed		= prio_changed_energy,
	.switched_to		= switched_to_energy,
};

__init void init_sched_energy_class(void)
{
	int i;
	hrtimer_init(&hr_timer, CLOCK_MONOTONIC, HRTIMER_MODE_REL);
	hr_timer.function = sched_period_timer;
	for (i = 0; i < _MAX_CPU; i++) 
		rq_lock[i] = __SPIN_LOCK_UNLOCKED(rq_lock[i]);
}
