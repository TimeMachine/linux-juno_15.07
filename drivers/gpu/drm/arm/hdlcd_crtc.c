/*
 * Copyright (C) 2013,2014 ARM Limited
 * Author: Liviu Dudau <Liviu.Dudau@arm.com>
 *
 * This file is subject to the terms and conditions of the GNU General Public
 * License.  See the file COPYING in the main directory of this archive
 * for more details.
 *
 *  Implementation of a CRTC class for the HDLCD driver.
 */

#include <linux/clk.h>
#include <drm/drmP.h>
#include <drm/drm_crtc.h>
#include <drm/drm_crtc_helper.h>
#include <drm/drm_fb_helper.h>
#include <drm/drm_fb_cma_helper.h>
#include <drm/drm_gem_cma_helper.h>

#include "hdlcd_drv.h"
#include "hdlcd_regs.h"

/*
 * The HDLCD controller is a dumb RGB streamer that gets connected to
 * a single HDMI transmitter or in the case of the ARM Models it gets
 * emulated by the software that does the actual rendering.
 *
 */
static void hdlcd_crtc_destroy(struct drm_crtc *crtc)
{
	drm_crtc_cleanup(crtc);
}

void hdlcd_set_scanout(struct hdlcd_drm_private *hdlcd, bool wait)
{
	struct drm_framebuffer *fb = hdlcd->crtc.fb;
	struct hdlcd_bo *bo;
	unsigned int depth, bpp;
	dma_addr_t scanout_start;
	int ret;

	drm_fb_get_bpp_depth(fb->pixel_format, &depth, &bpp);
	bo = hdlcd->bo;

	scanout_start = bo->dma_addr + fb->offsets[0] +
		(hdlcd->crtc.y * fb->pitches[0]) + (hdlcd->crtc.x * bpp/8);

	hdlcd_write(hdlcd, HDLCD_REG_FB_BASE, scanout_start);

	if (wait && hdlcd->dpms == DRM_MODE_DPMS_ON) {
		drm_vblank_get(fb->dev, 0);
		hdlcd->frame_completion.done = 0;
		do {
			ret = wait_for_completion_interruptible_timeout(&hdlcd->frame_completion,
							msecs_to_jiffies(50));
		} while (ret <= 0);
		drm_vblank_put(fb->dev, 0);
	} else {
		dev_info(fb->dev->dev, "%s: wait called with DPMS set to %d\n",
			__func__, hdlcd->dpms);
	}
}

static int hdlcd_crtc_page_flip(struct drm_crtc *crtc,
			struct drm_framebuffer *fb,
			struct drm_pending_vblank_event *event)
{
	struct hdlcd_drm_private *hdlcd = crtc_to_hdlcd_priv(crtc);

	if (hdlcd->dpms == DRM_MODE_DPMS_ON) {
		/* don't schedule any page flipping if one is in progress */
		if (hdlcd->event)
			return -EBUSY;

		hdlcd->event = event;
		drm_vblank_get(crtc->dev, 0);
	}

	crtc->fb = fb;

	if (hdlcd->dpms == DRM_MODE_DPMS_ON) {
		hdlcd_set_scanout(hdlcd, true);
	} else {
		unsigned long flags;

		/* not active, update registers immediately */
		hdlcd_set_scanout(hdlcd, false);
		spin_lock_irqsave(&crtc->dev->event_lock, flags);
		if (event)
			drm_send_vblank_event(crtc->dev, 0, event);
		spin_unlock_irqrestore(&crtc->dev->event_lock, flags);
	}

	return 0;
}

static const struct drm_crtc_funcs hdlcd_crtc_funcs = {
	.destroy	= hdlcd_crtc_destroy,
	.set_config	= drm_crtc_helper_set_config,
	.page_flip	= hdlcd_crtc_page_flip,
};

static void hdlcd_crtc_dpms(struct drm_crtc *crtc, int mode)
{
	struct hdlcd_drm_private *hdlcd = crtc_to_hdlcd_priv(crtc);

	hdlcd->dpms = mode;
	if (mode == DRM_MODE_DPMS_ON)
		hdlcd_write(hdlcd, HDLCD_REG_COMMAND, 1);
	else
		hdlcd_write(hdlcd, HDLCD_REG_COMMAND, 0);
}

static bool hdlcd_crtc_mode_fixup(struct drm_crtc *crtc,
			const struct drm_display_mode *mode,
			struct drm_display_mode *adjusted_mode)
{
	return true;
}

static void hdlcd_crtc_prepare(struct drm_crtc *crtc)
{
	drm_vblank_pre_modeset(crtc->dev, 0);
	hdlcd_crtc_dpms(crtc, DRM_MODE_DPMS_OFF);
}

static void hdlcd_crtc_commit(struct drm_crtc *crtc)
{
	drm_vblank_post_modeset(crtc->dev, 0);
	hdlcd_crtc_dpms(crtc, DRM_MODE_DPMS_ON);
}

static bool hdlcd_fb_mode_equal(struct drm_framebuffer *oldfb,
				struct drm_framebuffer *newfb)
{
	if (!oldfb || !newfb)
		return false;

	if (oldfb->pixel_format == newfb->pixel_format &&
		oldfb->width == newfb->width &&
		oldfb->height == newfb->height)
		return true;

	return false;
}

static int hdlcd_crtc_mode_set(struct drm_crtc *crtc,
			struct drm_display_mode *mode,
			struct drm_display_mode *adjusted_mode,
			int x, int y, struct drm_framebuffer *oldfb)
{
	struct hdlcd_drm_private *hdlcd = crtc_to_hdlcd_priv(crtc);
	unsigned int depth, bpp, polarities;
	unsigned char red_width = 0, green_width = 0, blue_width = 0, alpha_width = 0;
	unsigned int default_color = 0x00000000;

#ifdef HDLCD_SHOW_UNDERRUN
	default_color = 0x00ff000000;
#endif

	/* This function gets called when the only change is the start of
	   the scanout buffer. Detect that and bail out early */
	if (hdlcd->initialised && hdlcd_fb_mode_equal(oldfb, crtc->fb)) {
		hdlcd_set_scanout(hdlcd, true);
		return 0;
	}

	/* Preset the number of bits per colour */
	drm_fb_get_bpp_depth(crtc->fb->pixel_format, &depth, &bpp);
	switch (depth) {
	case 32:
		alpha_width = 8;
	case 24:
	case 8:	 /* pseudocolor */
		red_width = 8; green_width = 8; blue_width = 8;
		break;
	case 16: /* 565 format */
		red_width = 5; green_width = 6; blue_width = 5;
		break;
	}

	/* switch to using the more useful bytes per pixel */
	bpp = (bpp + 7) / 8;

	polarities = HDLCD_POLARITY_DATAEN | HDLCD_POLARITY_DATA;

	if (adjusted_mode->flags & DRM_MODE_FLAG_PHSYNC)
		polarities |= HDLCD_POLARITY_HSYNC;
	if (adjusted_mode->flags & DRM_MODE_FLAG_PVSYNC)
		polarities |= HDLCD_POLARITY_VSYNC;

	/* Allow max number of outstanding requests and largest burst size */
	hdlcd_write(hdlcd, HDLCD_REG_BUS_OPTIONS,
		    HDLCD_BUS_MAX_OUTSTAND | HDLCD_BUS_BURST_16);

	hdlcd_write(hdlcd, HDLCD_REG_PIXEL_FORMAT, (bpp - 1) << 3);

	hdlcd_write(hdlcd, HDLCD_REG_FB_LINE_LENGTH, crtc->fb->width * bpp);
	hdlcd_write(hdlcd, HDLCD_REG_FB_LINE_COUNT, crtc->fb->height - 1);
	hdlcd_write(hdlcd, HDLCD_REG_FB_LINE_PITCH, crtc->fb->width * bpp);
	hdlcd_write(hdlcd, HDLCD_REG_V_BACK_PORCH,
				mode->vtotal - mode->vsync_end - 1);
	hdlcd_write(hdlcd, HDLCD_REG_V_FRONT_PORCH,
				mode->vsync_start - mode->vdisplay - 1);
	hdlcd_write(hdlcd, HDLCD_REG_V_SYNC,
				mode->vsync_end - mode->vsync_start - 1);
	hdlcd_write(hdlcd, HDLCD_REG_V_DATA, mode->vdisplay - 1);
	hdlcd_write(hdlcd, HDLCD_REG_H_BACK_PORCH,
				mode->htotal - mode->hsync_end - 1);
	hdlcd_write(hdlcd, HDLCD_REG_H_FRONT_PORCH,
				mode->hsync_start - mode->hdisplay - 1);
	hdlcd_write(hdlcd, HDLCD_REG_H_SYNC,
				mode->hsync_end - mode->hsync_start - 1);
	hdlcd_write(hdlcd, HDLCD_REG_H_DATA, mode->hdisplay - 1);
	hdlcd_write(hdlcd, HDLCD_REG_POLARITIES, polarities);

	/*
	 * The format of the HDLCD_REG_<color>_SELECT register is:
	 *   - bits[23:16] - default value for that color component
	 *   - bits[11:8]  - number of bits to extract for each color component
	 *   - bits[4:0]   - index of the lowest bit to extract
	 *
	 * The default color value is used when bits[11:8] read zero, when the
	 * pixel is outside the visible frame area or when there is a
	 * buffer underrun.
	 */
	hdlcd_write(hdlcd, HDLCD_REG_BLUE_SELECT, default_color |
		alpha_width |   /* offset */
		(blue_width & 0xf) << 8);
	hdlcd_write(hdlcd, HDLCD_REG_GREEN_SELECT, default_color |
		(blue_width + alpha_width) |  /* offset */
		((green_width & 0xf) << 8));
	hdlcd_write(hdlcd, HDLCD_REG_RED_SELECT, default_color |
		(blue_width + green_width + alpha_width) |  /* offset */
		((red_width & 0xf) << 8));

	clk_prepare(hdlcd->clk);
	clk_set_rate(hdlcd->clk, mode->clock * 1000);
	clk_enable(hdlcd->clk);

	hdlcd_set_scanout(hdlcd, false);
	hdlcd->initialised = true;

	return 0;
}

int hdlcd_crtc_mode_set_base(struct drm_crtc *crtc, int x, int y,
			struct drm_framebuffer *oldfb)
{
	struct hdlcd_drm_private *hdlcd = crtc_to_hdlcd_priv(crtc);

	hdlcd_set_scanout(hdlcd, true);
	return 0;
}

static void hdlcd_crtc_load_lut(struct drm_crtc *crtc)
{
}

static const struct drm_crtc_helper_funcs hdlcd_crtc_helper_funcs = {
	.dpms		= hdlcd_crtc_dpms,
	.mode_fixup	= hdlcd_crtc_mode_fixup,
	.prepare	= hdlcd_crtc_prepare,
	.commit		= hdlcd_crtc_commit,
	.mode_set	= hdlcd_crtc_mode_set,
	.mode_set_base	= hdlcd_crtc_mode_set_base,
	.load_lut	= hdlcd_crtc_load_lut,
};

int hdlcd_setup_crtc(struct drm_device *dev)
{
	struct hdlcd_drm_private *hdlcd = dev->dev_private;
	int ret;

	drm_mode_config_init(dev);
	hdlcd_drm_mode_config_init(dev);

	ret = drm_crtc_init(dev, &hdlcd->crtc, &hdlcd_crtc_funcs);
	if (ret < 0)
		goto crtc_setup_err;

	drm_crtc_helper_add(&hdlcd->crtc, &hdlcd_crtc_helper_funcs);

	return 0;

crtc_setup_err:
	drm_mode_config_cleanup(dev);

	return ret;
}

