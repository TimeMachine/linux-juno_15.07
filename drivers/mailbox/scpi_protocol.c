/*
 * System Control and Power Interface (SCPI) Message Protocol driver
 *
 * SCPI Message Protocol is used between the System Control Processor(SCP)
 * and the Application Processors(AP). The Message Handling Unit(MHU)
 * provides a mechanism for inter-processor communication between SCP's
 * Cortex M3 and AP.
 *
 * SCP offers control and management of the core/cluster power states,
 * various power domain DVFS including the core/cluster, certain system
 * clocks configuration, thermal sensors and many others.
 *
 * Copyright (C) 2015 ARM Ltd.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms and conditions of the GNU General Public License,
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#define pr_fmt(fmt) KBUILD_MODNAME ": " fmt

#include <linux/bitmap.h>
#include <linux/device.h>
#include <linux/err.h>
#include <linux/export.h>
#include <linux/io.h>
#include <linux/kernel.h>
#include <linux/list.h>
#include <linux/mailbox_client.h>
#include <linux/module.h>
#include <linux/of_address.h>
#include <linux/of_platform.h>
#include <linux/printk.h>
#include <linux/scpi_protocol.h>
#include <linux/slab.h>
#include <linux/spinlock.h>

#define CMD_ID_SHIFT		0
#define CMD_ID_MASK		0x7f
#define CMD_TOKEN_ID_SHIFT	8
#define CMD_TOKEN_ID_MASK	0xff
#define CMD_DATA_SIZE_SHIFT	16
#define CMD_DATA_SIZE_MASK	0x1ff
#define PACK_SCPI_CMD(cmd_id, token, tx_sz)			\
	((((cmd_id) & CMD_ID_MASK) << CMD_ID_SHIFT) |		\
	(((token) & CMD_TOKEN_ID_MASK) << CMD_TOKEN_ID_SHIFT) |	\
	(((tx_sz) & CMD_DATA_SIZE_MASK) << CMD_DATA_SIZE_SHIFT))

#define CMD_SIZE(cmd)	(((cmd) >> CMD_DATA_SIZE_SHIFT) & CMD_DATA_SIZE_MASK)
#define CMD_UNIQ_MASK	(CMD_TOKEN_ID_MASK << CMD_TOKEN_ID_SHIFT | CMD_ID_MASK)
#define CMD_XTRACT_UNIQ(cmd)	((cmd) & CMD_UNIQ_MASK)

#define SCPI_SLOT		0

#define MAX_DVFS_DOMAINS	3
#define MAX_DVFS_OPPS		8
#define DVFS_LATENCY(hdr)	(le32_to_cpu(hdr) >> 16)
#define DVFS_OPP_COUNT(hdr)	((le32_to_cpu(hdr) >> 8) & 0xff)

#define PROTOCOL_REV_MINOR_BITS	16
#define PROTOCOL_REV_MINOR_MASK	((1U << PROTOCOL_REV_MINOR_BITS) - 1)
#define PROTOCOL_REV_MAJOR(x)	((x) >> PROTOCOL_REV_MINOR_BITS)
#define PROTOCOL_REV_MINOR(x)	((x) & PROTOCOL_REV_MINOR_MASK)

#define FW_REV_MAJOR_BITS	24
#define FW_REV_MINOR_BITS	16
#define FW_REV_PATCH_MASK	((1U << FW_REV_MINOR_BITS) - 1)
#define FW_REV_MINOR_MASK	((1U << FW_REV_MAJOR_BITS) - 1)
#define FW_REV_MAJOR(x)		((x) >> FW_REV_MAJOR_BITS)
#define FW_REV_MINOR(x)		(((x) & FW_REV_MINOR_MASK) >> FW_REV_MINOR_BITS)
#define FW_REV_PATCH(x)		((x) & FW_REV_PATCH_MASK)

#define MAX_RX_TIMEOUT		(msecs_to_jiffies(30))

enum scpi_error_codes {
	SCPI_SUCCESS = 0, /* Success */
	SCPI_ERR_PARAM = 1, /* Invalid parameter(s) */
	SCPI_ERR_ALIGN = 2, /* Invalid alignment */
	SCPI_ERR_SIZE = 3, /* Invalid size */
	SCPI_ERR_HANDLER = 4, /* Invalid handler/callback */
	SCPI_ERR_ACCESS = 5, /* Invalid access/permission denied */
	SCPI_ERR_RANGE = 6, /* Value out of range */
	SCPI_ERR_TIMEOUT = 7, /* Timeout has occurred */
	SCPI_ERR_NOMEM = 8, /* Invalid memory area or pointer */
	SCPI_ERR_PWRSTATE = 9, /* Invalid power state */
	SCPI_ERR_SUPPORT = 10, /* Not supported or disabled */
	SCPI_ERR_DEVICE = 11, /* Device error */
	SCPI_ERR_MAX
};

enum scpi_std_cmd {
	SCPI_CMD_INVALID		= 0x00,
	SCPI_CMD_SCPI_READY		= 0x01,
	SCPI_CMD_SCPI_CAPABILITIES	= 0x02,
	SCPI_CMD_SET_CSS_PWR_STATE	= 0x03,
	SCPI_CMD_GET_CSS_PWR_STATE	= 0x04,
	SCPI_CMD_SET_SYS_PWR_STATE	= 0x05,
	SCPI_CMD_SET_CPU_TIMER		= 0x06,
	SCPI_CMD_CANCEL_CPU_TIMER	= 0x07,
	SCPI_CMD_DVFS_CAPABILITIES	= 0x08,
	SCPI_CMD_GET_DVFS_INFO		= 0x09,
	SCPI_CMD_SET_DVFS		= 0x0a,
	SCPI_CMD_GET_DVFS		= 0x0b,
	SCPI_CMD_GET_DVFS_STAT		= 0x0c,
	SCPI_CMD_CLOCK_CAPABILITIES	= 0x0d,
	SCPI_CMD_GET_CLOCK_INFO		= 0x0e,
	SCPI_CMD_SET_CLOCK_VALUE	= 0x0f,
	SCPI_CMD_GET_CLOCK_VALUE	= 0x10,
	SCPI_CMD_PSU_CAPABILITIES	= 0x11,
	SCPI_CMD_GET_PSU_INFO		= 0x12,
	SCPI_CMD_SET_PSU		= 0x13,
	SCPI_CMD_GET_PSU		= 0x14,
	SCPI_CMD_SENSOR_CAPABILITIES	= 0x15,
	SCPI_CMD_SENSOR_INFO		= 0x16,
	SCPI_CMD_SENSOR_VALUE		= 0x17,
	SCPI_CMD_SENSOR_CFG_PERIODIC	= 0x18,
	SCPI_CMD_SENSOR_CFG_BOUNDS	= 0x19,
	SCPI_CMD_SENSOR_ASYNC_VALUE	= 0x1a,
	SCPI_CMD_SET_DEVICE_PWR_STATE	= 0x1b,
	SCPI_CMD_GET_DEVICE_PWR_STATE	= 0x1c,
	SCPI_CMD_COUNT
};

struct scpi_xfer {
	u32 slot; /* has to be first element */
	u32 cmd;
	u32 status;
	const void *tx_buf;
	void *rx_buf;
	unsigned int tx_len;
	struct list_head node;
	struct completion done;
};

struct scpi_chan {
	struct mbox_client cl;
	struct mbox_chan *chan;
	void __iomem *tx_payload;
	void __iomem *rx_payload;
	struct list_head rx_pending;
	struct list_head xfers_list;
	struct scpi_xfer *xfers;
	spinlock_t rx_lock; /* locking for the rx pending list */
	struct mutex xfers_lock;
	atomic_t token;
};

struct scpi_drvinfo {
	u32 protocol_version;
	u32 firmware_version;
	int num_chans;
	atomic_t next_chan;
	struct scpi_ops *scpi_ops;
	struct scpi_chan *channels;
	struct scpi_dvfs_info *dvfs[MAX_DVFS_DOMAINS];
};

/*
 * The SCP firmware only executes in little-endian mode, so any buffers
 * shared through SCPI should have their contents converted to little-endian
 */
struct scpi_shared_mem {
	__le32 command;
	__le32 status;
	u8 payload[0];
} __packed;

struct scp_capabilities {
	__le32 protocol_version;
	__le32 event_version;
	__le32 platform_version;
	__le32 commands[4];
} __packed;

struct clk_get_info {
	__le16 id;
	__le16 flags;
	__le32 min_rate;
	__le32 max_rate;
	u8 name[20];
} __packed;

struct clk_get_value {
	__le32 rate;
} __packed;

struct clk_set_value {
	__le16 id;
	__le16 reserved;
	__le32 rate;
} __packed;

struct dvfs_info {
	__le32 header;
	struct {
		__le32 freq;
		__le32 m_volt;
	} opps[MAX_DVFS_OPPS];
} __packed;

struct dvfs_get {
	u8 index;
} __packed;

struct dvfs_set {
	u8 domain;
	u8 index;
} __packed;

struct sensor_capabilities {
	__le16 sensors;
} __packed;

struct sensor_info {
	__le16 id;
	u8 class;
	u8 triggers;
	u8 name[20];
} __packed;

struct sensor_value {
	__le32 value;
} __packed;

static struct scpi_drvinfo *scpi_info;

static int scpi_linux_errmap[SCPI_ERR_MAX] = {
	/* better than switch case as long as return value is continuous */
	0, /* SCPI_SUCCESS */
	-EINVAL, /* SCPI_ERR_PARAM */
	-ENOEXEC, /* SCPI_ERR_ALIGN */
	-EMSGSIZE, /* SCPI_ERR_SIZE */
	-EINVAL, /* SCPI_ERR_HANDLER */
	-EACCES, /* SCPI_ERR_ACCESS */
	-ERANGE, /* SCPI_ERR_RANGE */
	-ETIMEDOUT, /* SCPI_ERR_TIMEOUT */
	-ENOMEM, /* SCPI_ERR_NOMEM */
	-EINVAL, /* SCPI_ERR_PWRSTATE */
	-EOPNOTSUPP, /* SCPI_ERR_SUPPORT */
	-EIO, /* SCPI_ERR_DEVICE */
};

static inline int scpi_to_linux_errno(int errno)
{
	if (errno >= SCPI_SUCCESS && errno < SCPI_ERR_MAX)
		return scpi_linux_errmap[errno];
	return -EIO;
}

/**
 * sysfs_create_groups - given a directory kobject, create a bunch of attribute groups
 * @kobj:	The kobject to create the group on
 * @groups:	The attribute groups to create, NULL terminated
 *
 * This function creates a bunch of attribute groups.  If an error occurs when
 * creating a group, all previously created groups will be removed, unwinding
 * everything back to the original state when this function was called.
 * It will explicitly warn and error if any of the attribute files being
 * created already exist.
 *
 * Returns 0 on success or error code from sysfs_create_groups on error.
 */
static int sysfs_create_groups(struct kobject *kobj,
			const struct attribute_group **groups)
{
	int error = 0;
	int i;

	if (!groups)
		return 0;

	for (i = 0; groups[i]; i++) {
		error = sysfs_create_group(kobj, groups[i]);
		if (error) {
			while (--i >= 0)
				sysfs_remove_group(kobj, groups[i]);
			break;
		}
	}
	return error;
}

/**
 * sysfs_remove_groups - remove a list of groups
 *
 * kobj:	The kobject for the groups to be removed from
 * groups:	NULL terminated list of groups to be removed
 *
 * If groups is not NULL, the all groups will be removed from the kobject
 */
static void sysfs_remove_groups(struct kobject *kobj,
			 const struct attribute_group **groups)
{
	int i;

	if (!groups)
		return;
	for (i = 0; groups[i]; i++)
		sysfs_remove_group(kobj, groups[i]);
}

static void scpi_process_cmd(struct scpi_chan *ch, u32 cmd)
{
	unsigned long flags;
	struct scpi_xfer *t, *match = NULL;

	spin_lock_irqsave(&ch->rx_lock, flags);
	if (list_empty(&ch->rx_pending)) {
		spin_unlock_irqrestore(&ch->rx_lock, flags);
		return;
	}

	list_for_each_entry(t, &ch->rx_pending, node)
		if (CMD_XTRACT_UNIQ(t->cmd) == CMD_XTRACT_UNIQ(cmd)) {
			list_del(&t->node);
			match = t;
			break;
		}
	/* check if wait_for_completion is in progress or timed-out */
	if (match && !completion_done(&match->done)) {
		struct scpi_shared_mem *mem = ch->rx_payload;

		match->status = le32_to_cpu(mem->status);
		memcpy_fromio(match->rx_buf, mem->payload, CMD_SIZE(cmd));
		complete(&match->done);
	}
	spin_unlock_irqrestore(&ch->rx_lock, flags);
}

static void scpi_handle_remote_msg(struct mbox_client *c, void *msg)
{
	struct scpi_chan *ch = container_of(c, struct scpi_chan, cl);
	struct scpi_shared_mem *mem = ch->rx_payload;
	u32 cmd = le32_to_cpu(mem->command);

	scpi_process_cmd(ch, cmd);
}

static void scpi_tx_prepare(struct mbox_client *c, void *msg)
{
	unsigned long flags;
	struct scpi_xfer *t = msg;
	struct scpi_chan *ch = container_of(c, struct scpi_chan, cl);
	struct scpi_shared_mem *mem = (struct scpi_shared_mem *)ch->tx_payload;

	mem->command = cpu_to_le32(t->cmd);
	if (t->tx_buf)
		memcpy_toio(mem->payload, t->tx_buf, t->tx_len);
	if (t->rx_buf) {
		spin_lock_irqsave(&ch->rx_lock, flags);
		list_add_tail(&t->node, &ch->rx_pending);
		spin_unlock_irqrestore(&ch->rx_lock, flags);
	}
}

static struct scpi_xfer *get_scpi_xfer(struct scpi_chan *ch)
{
	struct scpi_xfer *t;

	mutex_lock(&ch->xfers_lock);
	if (list_empty(&ch->xfers_list)) {
		mutex_unlock(&ch->xfers_lock);
		return NULL;
	}
	t = list_first_entry(&ch->xfers_list, struct scpi_xfer, node);
	list_del(&t->node);
	mutex_unlock(&ch->xfers_lock);
	return t;
}

static void put_scpi_xfer(struct scpi_xfer *t, struct scpi_chan *ch)
{
	mutex_lock(&ch->xfers_lock);
	list_add_tail(&t->node, &ch->xfers_list);
	mutex_unlock(&ch->xfers_lock);
}

static int
scpi_send_message(u8 cmd, void *tx_buf, unsigned int len, void *rx_buf)
{
	int ret;
	u8 token, chan;
	struct scpi_xfer *msg;
	struct scpi_chan *scpi_chan;

	chan = atomic_inc_return(&scpi_info->next_chan) % scpi_info->num_chans;
	scpi_chan = scpi_info->channels + chan;

	msg = get_scpi_xfer(scpi_chan);
	if (!msg)
		return -ENOMEM;

	token = atomic_inc_return(&scpi_chan->token) & CMD_TOKEN_ID_MASK;

	msg->slot = BIT(SCPI_SLOT);
	msg->cmd = PACK_SCPI_CMD(cmd, token, len);
	msg->tx_buf = tx_buf;
	msg->tx_len = len;
	msg->rx_buf = rx_buf;
	init_completion(&msg->done);

	ret = mbox_send_message(scpi_chan->chan, msg);
	if (ret < 0 || !rx_buf)
		goto out;

	if (!wait_for_completion_timeout(&msg->done, MAX_RX_TIMEOUT))
		ret = -ETIMEDOUT;
	else
		/* first status word */
		ret = le32_to_cpu(msg->status);
out:
	if (ret < 0 && rx_buf) /* remove entry from the list if timed-out */
		scpi_process_cmd(scpi_chan, msg->cmd);

	put_scpi_xfer(msg, scpi_chan);
	/* SCPI error codes > 0, translate them to Linux scale*/
	return ret > 0 ? scpi_to_linux_errno(ret) : ret;
}

static u32 scpi_get_version(void)
{
	return scpi_info->protocol_version;
}

static int
scpi_clk_get_range(u16 clk_id, unsigned long *min, unsigned long *max)
{
	int ret;
	struct clk_get_info clk;
	__le16 le_clk_id = cpu_to_le16(clk_id);

	ret = scpi_send_message(SCPI_CMD_GET_CLOCK_INFO,
				&le_clk_id, sizeof(le_clk_id), &clk);
	if (!ret) {
		*min = le32_to_cpu(clk.min_rate);
		*max = le32_to_cpu(clk.max_rate);
	}
	return ret;
}

static unsigned long scpi_clk_get_val(u16 clk_id)
{
	int ret;
	struct clk_get_value clk;
	__le16 le_clk_id = cpu_to_le16(clk_id);

	ret = scpi_send_message(SCPI_CMD_GET_CLOCK_VALUE,
				&le_clk_id, sizeof(le_clk_id), &clk);
	return ret ? ret : le32_to_cpu(clk.rate);
}

static int scpi_clk_set_val(u16 clk_id, unsigned long rate)
{
	int stat;
	struct clk_set_value clk = { cpu_to_le16(clk_id), 0, cpu_to_le32(rate) };

	return scpi_send_message(SCPI_CMD_SET_CLOCK_VALUE,
				 &clk, sizeof(clk), &stat);
}

static int scpi_dvfs_get_idx(u8 domain)
{
	int ret;
	struct dvfs_get dvfs;

	ret = scpi_send_message(SCPI_CMD_GET_DVFS,
				&domain, sizeof(domain), &dvfs);
	return ret ? ret : dvfs.index;
}

static int scpi_dvfs_set_idx(u8 domain, u8 index)
{
	int stat;
	struct dvfs_set dvfs = {domain, index};

	return scpi_send_message(SCPI_CMD_SET_DVFS,
				 &dvfs, sizeof(dvfs), &stat);
}

static struct scpi_dvfs_info *scpi_dvfs_get_info(u8 domain)
{
	struct scpi_dvfs_info *info;
	struct scpi_opp *opp;
	struct dvfs_info buf;
	int ret, i;

	if (domain >= MAX_DVFS_DOMAINS)
		return ERR_PTR(-EINVAL);

	if (scpi_info->dvfs[domain])	/* data already populated */
		return scpi_info->dvfs[domain];

	ret = scpi_send_message(SCPI_CMD_GET_DVFS_INFO, &domain,
				sizeof(domain), &buf);

	if (ret)
		return ERR_PTR(ret);

	info = kmalloc(sizeof(*info), GFP_KERNEL);
	if (!info)
		return ERR_PTR(-ENOMEM);

	info->count = DVFS_OPP_COUNT(buf.header);
	info->latency = DVFS_LATENCY(buf.header) * 1000; /* uS to nS */

	info->opps = kcalloc(info->count, sizeof(*opp), GFP_KERNEL);
	if (!info->opps) {
		kfree(info);
		return ERR_PTR(-ENOMEM);
	}

	for (i = 0, opp = info->opps; i < info->count; i++, opp++) {
		opp->freq = le32_to_cpu(buf.opps[i].freq);
		opp->m_volt = le32_to_cpu(buf.opps[i].m_volt);
	}

	scpi_info->dvfs[domain] = info;
	return info;
}

static int scpi_sensor_get_id(char *name)
{
	struct sensor_capabilities caps;
	u16 sensors, sensor_id;
	int ret;

	ret = scpi_send_message(SCPI_CMD_SENSOR_CAPABILITIES,
				NULL, 0, &caps);
	if (ret)
		return ret;
	sensors = le16_to_cpu(caps.sensors);

	for (sensor_id = 0; sensor_id < sensors; sensor_id++) {

		struct sensor_info info;
		__le16 le_sensor_id = cpu_to_le16(sensor_id);
		ret = scpi_send_message(SCPI_CMD_SENSOR_INFO, &le_sensor_id,
					sizeof(le_sensor_id), &info);
		if (ret)
			return ret;

		if (strcmp(name, info.name) == 0)
			return sensor_id;
	}

	return -ENODEV;
}

static int scpi_sensor_get_value(u16 sensor_id, u32 *value)
{
	struct sensor_value val;
	__le16 le_sensor_id = cpu_to_le16(sensor_id);
	int ret;

	ret = scpi_send_message(SCPI_CMD_SENSOR_VALUE, &le_sensor_id,
				sizeof(le_sensor_id), &val);
	if (ret == 0)
		*value = le32_to_cpu(val.value);

	return ret;
}

static struct scpi_ops scpi_ops = {
	.get_version = scpi_get_version,
	.clk_get_range = scpi_clk_get_range,
	.clk_get_val = scpi_clk_get_val,
	.clk_set_val = scpi_clk_set_val,
	.dvfs_get_idx = scpi_dvfs_get_idx,
	.dvfs_set_idx = scpi_dvfs_set_idx,
	.dvfs_get_info = scpi_dvfs_get_info,
	.sensor_get_id = scpi_sensor_get_id,
	.sensor_get_value = scpi_sensor_get_value,
};

struct scpi_ops *get_scpi_ops(void)
{
	return scpi_info ? scpi_info->scpi_ops : NULL;
}
EXPORT_SYMBOL_GPL(get_scpi_ops);

static int scpi_init_versions(struct scpi_drvinfo *info)
{
	int ret;
	struct scp_capabilities caps;

	ret = scpi_send_message(SCPI_CMD_SCPI_CAPABILITIES, NULL, 0, &caps);
	if (!ret) {
		info->protocol_version = le32_to_cpu(caps.protocol_version);
		info->firmware_version = le32_to_cpu(caps.platform_version);
	}
	return ret;
}

static ssize_t protocol_version_show(struct device *dev,
				     struct device_attribute *attr, char *buf)
{
	struct scpi_drvinfo *scpi_info = dev_get_drvdata(dev);

	return sprintf(buf, "%d.%d\n",
		       PROTOCOL_REV_MAJOR(scpi_info->protocol_version),
		       PROTOCOL_REV_MINOR(scpi_info->protocol_version));
}
static DEVICE_ATTR_RO(protocol_version);

static ssize_t firmware_version_show(struct device *dev,
				     struct device_attribute *attr, char *buf)
{
	struct scpi_drvinfo *scpi_info = dev_get_drvdata(dev);

	return sprintf(buf, "%d.%d.%d\n",
		       FW_REV_MAJOR(scpi_info->firmware_version),
		       FW_REV_MINOR(scpi_info->firmware_version),
		       FW_REV_PATCH(scpi_info->firmware_version));
}
static DEVICE_ATTR_RO(firmware_version);

static struct attribute *versions_attrs[] = {
	&dev_attr_firmware_version.attr,
	&dev_attr_protocol_version.attr,
	NULL,
};
ATTRIBUTE_GROUPS(versions);

static void
scpi_free_channels(struct device *dev, struct scpi_chan *pchan, int count)
{
	int i;

	for (i = 0; i < count && pchan->chan; i++, pchan++) {
		mbox_free_channel(pchan->chan);
		devm_kfree(dev, pchan->xfers);
		devm_iounmap(dev, pchan->rx_payload);
	}
}

static int scpi_remove(struct platform_device *pdev)
{
	struct device *dev = &pdev->dev;
	struct scpi_drvinfo *scpi_info = platform_get_drvdata(pdev);

	of_platform_depopulate(dev);
	sysfs_remove_groups(&dev->kobj, versions_groups);
	scpi_free_channels(dev, scpi_info->channels, scpi_info->num_chans);
	platform_set_drvdata(pdev, NULL);

	devm_kfree(dev, scpi_info->channels);
	devm_kfree(dev, scpi_info);
	scpi_info = NULL;

	return 0;
}

#define MAX_SCPI_XFERS		10
static int scpi_alloc_xfer_list(struct device *dev, struct scpi_chan *ch)
{
	int i;
	struct scpi_xfer *xfers;

	xfers = devm_kzalloc(dev, MAX_SCPI_XFERS * sizeof(*xfers), GFP_KERNEL);
	if (!xfers)
		return -ENOMEM;

	ch->xfers = xfers;
	for (i = 0; i < MAX_SCPI_XFERS; i++, xfers++)
		list_add_tail(&xfers->node, &ch->xfers_list);
	return 0;
}

static int scpi_probe(struct platform_device *pdev)
{
	int count, idx, ret;
	struct resource res;
	struct scpi_chan *scpi_chan;
	struct device *dev = &pdev->dev;
	struct device_node *np = dev->of_node;

	scpi_info = devm_kzalloc(dev, sizeof(*scpi_info), GFP_KERNEL);
	if (!scpi_info) {
		dev_err(dev, "failed to allocate memory for scpi drvinfo\n");
		return -ENOMEM;
	}

	count = of_count_phandle_with_args(np, "mboxes", "#mbox-cells");
	if (count < 0) {
		dev_err(dev, "no mboxes property in '%s'\n", np->full_name);
		return -ENODEV;
	}

	scpi_chan = devm_kcalloc(dev, count, sizeof(*scpi_chan), GFP_KERNEL);
	if (!scpi_chan) {
		dev_err(dev, "failed to allocate memory scpi chaninfo\n");
		return -ENOMEM;
	}

	for (idx = 0; idx < count; idx++) {
		resource_size_t size;
		struct scpi_chan *pchan = scpi_chan + idx;
		struct mbox_client *cl = &pchan->cl;
		struct device_node *shmem = of_parse_phandle(np, "shmem", idx);

		if (of_address_to_resource(shmem, 0, &res)) {
			dev_err(dev, "failed to get SCPI payload mem resource\n");
			ret = -EINVAL;
			goto err;
		}

		size = resource_size(&res);
		pchan->rx_payload = devm_ioremap(dev, res.start, size);
		if (!pchan->rx_payload) {
			dev_err(dev, "failed to ioremap SCPI payload\n");
			ret = -EADDRNOTAVAIL;
			goto err;
		}
		pchan->tx_payload = pchan->rx_payload + (size >> 1);

		cl->dev = dev;
		cl->rx_callback = scpi_handle_remote_msg;
		cl->tx_prepare = scpi_tx_prepare;
		cl->tx_block = true;
		cl->tx_tout = 50;
		cl->knows_txdone = false; /* controller can ack */

		INIT_LIST_HEAD(&pchan->rx_pending);
		INIT_LIST_HEAD(&pchan->xfers_list);
		spin_lock_init(&pchan->rx_lock);
		mutex_init(&pchan->xfers_lock);

		ret = scpi_alloc_xfer_list(dev, pchan);
		if (!ret) {
			pchan->chan = mbox_request_channel(cl, idx);
			if (!IS_ERR(pchan->chan))
				continue;
			ret = -EPROBE_DEFER;
			dev_err(dev, "failed to acquire channel#%d\n", idx);
		}
err:
		scpi_free_channels(dev, scpi_chan, idx);
		scpi_info = NULL;
		return ret;
	}

	scpi_info->channels = scpi_chan;
	scpi_info->num_chans = count;
	platform_set_drvdata(pdev, scpi_info);

	ret = scpi_init_versions(scpi_info);
	if (ret) {
		dev_err(dev, "incorrect or no SCP firmware found\n");
		scpi_remove(pdev);
		return ret;
	}

	_dev_info(dev, "SCP Protocol %d.%d Firmware %d.%d.%d version\n",
		  PROTOCOL_REV_MAJOR(scpi_info->protocol_version),
		  PROTOCOL_REV_MINOR(scpi_info->protocol_version),
		  FW_REV_MAJOR(scpi_info->firmware_version),
		  FW_REV_MINOR(scpi_info->firmware_version),
		  FW_REV_PATCH(scpi_info->firmware_version));
	scpi_info->scpi_ops = &scpi_ops;

	ret = sysfs_create_groups(&dev->kobj, versions_groups);
	if (ret)
		dev_err(dev, "unable to create sysfs version group\n");

	return of_platform_populate(dev->of_node, NULL, NULL, dev);
}

static const struct of_device_id scpi_of_match[] = {
	{.compatible = "arm,scpi"},
	{},
};

MODULE_DEVICE_TABLE(of, scpi_of_match);

static struct platform_driver scpi_driver = {
	.driver = {
		.name = "scpi_protocol",
		.of_match_table = scpi_of_match,
	},
	.probe = scpi_probe,
	.remove = scpi_remove,
};
module_platform_driver(scpi_driver);

MODULE_AUTHOR("Sudeep Holla <sudeep.holla@arm.com>");
MODULE_DESCRIPTION("ARM SCPI mailbox protocol driver");
MODULE_LICENSE("GPL");
