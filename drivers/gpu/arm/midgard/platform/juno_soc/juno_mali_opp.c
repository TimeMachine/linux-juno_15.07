/*
 *
 * (C) COPYRIGHT ARM Limited. All rights reserved.
 *
 * This program is free software and is provided to you under the terms of the
 * GNU General Public License version 2 as published by the Free Software
 * Foundation, and any use by you of this program is subject to the terms
 * of such GNU licence.
 *
 * A copy of the licence is included with the program, and can also be obtained
 * from Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 *
 */



#include <linux/module.h>
#include <linux/of_platform.h>
#include <linux/platform_device.h>
#include <linux/opp.h>
#include <linux/scpi_protocol.h>

#include <mali_kbase.h>

static int init_juno_opps_from_scpi(struct device *dev)
{
	struct scpi_ops *scpi;
	struct scpi_dvfs_info *info;
	int i;

	scpi = get_scpi_ops();
	if (!scpi)
		return 0; /* Really need to defer until scpi available */

	/* Hard coded for Juno. 2 is GPU domain */
	info = scpi->dvfs_get_info(2);
	if (IS_ERR_OR_NULL(info))
		return PTR_ERR(info);

	for (i = 0; i < info->count; i++) {
		struct scpi_opp *e = &info->opps[i];
		dev_info(dev, "Mali OPP from SCPI: %u Hz @ %u mV\n",
				e->freq, e->m_volt);
		opp_add(dev, e->freq, e->m_volt * 1000);
	}

	return 0;
}

int setup_opps(void)
{
	struct device_node *np;
	struct platform_device *pdev;
	int err;

	np = of_find_node_by_name(NULL, "gpu");
	if (!np) {
		printk(KERN_ERR "Failed to find DT entry for Mali\n");
		return -EFAULT;
	}

	pdev = of_find_device_by_node(np);
	if (!pdev) {
		printk(KERN_ERR "Failed to find device for Mali\n");
		of_node_put(np);
		return -EFAULT;
	}

	err = init_juno_opps_from_scpi(&pdev->dev);

	of_node_put(np);

	return err;
}
