#ifndef _MP4_OUTPUT_H
#define _MP4_OUTPUT_H

#include "core/siril.h"

#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
#include <libavutil/frame.h>

/* same order as in the combo box 'comboExport' */
typedef enum {
	EXPORT_FITS,
	EXPORT_FITSEQ,
	EXPORT_TIFF,
	EXPORT_SER,
	EXPORT_AVI,
	EXPORT_MP4,
	EXPORT_MP4_H265,
	EXPORT_WEBM
} export_format;

struct mp4_struct {
	AVOutputFormat *fmt;
	AVFormatContext *oc;
	AVStream *st;
	AVCodecContext *enc;

	/* pts of the next frame that will be generated */
	int64_t next_pts;
	int samples_count;

	AVFrame *frame;
	AVFrame *tmp_frame;

	struct SwsContext *sws_ctx;
	struct SwrContext *swr_ctx;

	int64_t bitrate;
	int src_w, src_h;

};

struct mp4_struct *mp4_create(const char *filename, int dst_w, int dst_h, int fps, int nb_layers, int quality, int src_w, int src_h, export_format type);
int mp4_add_frame(struct mp4_struct *, fits *);
int mp4_close(struct mp4_struct *);

#endif
