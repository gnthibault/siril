/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_FFMPEG

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <libavutil/avassert.h>
#include <libavutil/channel_layout.h>
#include <libavutil/opt.h>
#include <libavutil/mathematics.h>
#include <libavutil/timestamp.h>
#include <libswscale/swscale.h>
#include <libswresample/swresample.h>

#include "io/mp4_output.h"
#include "core/proto.h"	// computePente
#include "gui/callbacks.h"	// siril_log_message

//#define STREAM_DURATION   10.0
#define STREAM_PIX_FMT    AV_PIX_FMT_YUV420P /* default output pix_fmt */

#define SCALE_FLAGS SWS_BICUBIC

/* Formats and codecs:
 * mp4 output should use AV_CODEC_ID_H264 codec
 * webm output should use AV_CODEC_ID_VP8 codec, default is VP9
 */


static void log_packet(const AVFormatContext *fmt_ctx, const AVPacket *pkt)
{
	AVRational *time_base = &fmt_ctx->streams[pkt->stream_index]->time_base;

	printf("pts:%s pts_time:%s dts:%s dts_time:%s duration:%s duration_time:%s stream_index:%d\n",
			av_ts2str(pkt->pts), av_ts2timestr(pkt->pts, time_base),
			av_ts2str(pkt->dts), av_ts2timestr(pkt->dts, time_base),
			av_ts2str(pkt->duration), av_ts2timestr(pkt->duration, time_base),
			pkt->stream_index);
}

static int write_frame(AVFormatContext *fmt_ctx, const AVRational *time_base, AVStream *st, AVPacket *pkt)
{
	/* rescale output packet timestamp values from codec to stream timebase */
	av_packet_rescale_ts(pkt, *time_base, st->time_base);
	pkt->stream_index = st->index;

	/* Write the compressed frame to the media file. */
	log_packet(fmt_ctx, pkt);
	return av_interleaved_write_frame(fmt_ctx, pkt);
	// avcodec_receive_packet?
}

/* Add an output stream. */
static int add_stream(struct mp4_struct *ost, AVCodec **codec,
		enum AVCodecID codec_id, int w, int h, int fps)
{
	AVCodecContext *c;

	/* find the encoder */
	*codec = avcodec_find_encoder(codec_id);
	if (!(*codec)) {
		fprintf(stderr, "Could not find encoder for '%s'\n", avcodec_get_name(codec_id));
		return 1;
	}

	ost->st = avformat_new_stream(ost->oc, NULL);
	if (!ost->st) {
		fprintf(stderr, "Could not allocate stream\n");
		return 1;
	}
	ost->st->id = ost->oc->nb_streams-1;
	c = avcodec_alloc_context3(*codec);
	if (!c) {
		fprintf(stderr, "Could not alloc an encoding context\n");
		return 1;
	}
	ost->enc = c;

	switch ((*codec)->type) {
		case AVMEDIA_TYPE_VIDEO:
			c->codec_id = codec_id;

			c->bit_rate = ost->bitrate;
			c->bit_rate_tolerance = 50000;
			/* Resolution must be a multiple of two. */
			c->width    = w;
			c->height   = h;
			/* timebase: This is the fundamental unit of time (in seconds) in terms
			 * of which frame timestamps are represented. For fixed-fps content,
			 * timebase should be 1/framerate and timestamp increments should be
			 * identical to 1. */
			ost->st->time_base = (AVRational){ 1, fps };
			c->time_base       = ost->st->time_base;

			c->gop_size      = 12; /* emit one intra frame every twelve frames at most */
			c->pix_fmt       = STREAM_PIX_FMT;
			break;

		default:
			fprintf(stderr, "strange case happening in stream addition\n");
			break;
	}

	/* Some formats want stream headers to be separate. */
	if (ost->oc->oformat->flags & AVFMT_GLOBALHEADER)
		c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

	return 0;
}

/**************************************************************/
/* video output */

static AVFrame *alloc_picture(enum AVPixelFormat pix_fmt, int width, int height)
{
	AVFrame *picture;
	int ret;

	picture = av_frame_alloc();
	if (!picture)
		return NULL;

	picture->format = pix_fmt;
	picture->width  = width;
	picture->height = height;

	/* allocate the buffers for the frame data */
	ret = av_frame_get_buffer(picture, 32);
	if (ret < 0) {
		fprintf(stderr, "Could not allocate frame data.\n");
		return NULL;
	}

	return picture;
}

static int open_video(AVCodec *codec, struct mp4_struct *ost, AVDictionary *opt_arg, int nb_layers)
{
	int ret;
	AVCodecContext *c = ost->enc;
	AVDictionary *opt = NULL;

	av_dict_copy(&opt, opt_arg, 0);

	/* open the codec */
	ret = avcodec_open2(c, codec, &opt);
	av_dict_free(&opt);
	if (ret < 0) {
		fprintf(stderr, "Could not open video codec: %s\n", av_err2str(ret));
		return 1;
	}

	/* allocate and init a re-usable frame */
	ost->frame = alloc_picture(c->pix_fmt, c->width, c->height);
	if (!ost->frame) {
		fprintf(stderr, "Could not allocate video frame\n");
		return 1;
	}

	/* If the output format is not RGB24, then a temporary RGB24 picture is
	 * needed too. It is then converted to the required output format. */
	ost->tmp_frame = NULL;
	if (c->pix_fmt != AV_PIX_FMT_RGB24) {
		enum AVPixelFormat pix_fmt = (nb_layers == 1) ? AV_PIX_FMT_GRAY8 : AV_PIX_FMT_RGB24;
		ost->tmp_frame = alloc_picture(pix_fmt, ost->src_w, ost->src_h);
		if (!ost->tmp_frame) {
			fprintf(stderr, "Could not allocate temporary picture\n");
			return 1;
		}
	}

	/* copy the stream parameters to the muxer */
	ret = avcodec_parameters_from_context(ost->st->codecpar, c);
	if (ret < 0) {
		fprintf(stderr, "Could not copy the stream parameters\n");
		return 1;
	}

	return 0;
}

/* void RGB2YUV(WORD r, WORD g, WORD b, BYTE *y, BYTE *u, BYTE *v) {
	double R = (double)r, G = (double)g, B = (double)b;
	*y = round_to_WORD((0.257 * R) + (0.504 * G) + (0.098 * B) + 16.0);
	*v = round_to_WORD((0.439 * R) - (0.368 * G) - (0.071 * B) + 128.0);
	*u = round_to_WORD(-(0.148 * R) - (0.291 * G) + (0.439 * B) + 128.0);
}*/

static int fill_rgb_image(AVFrame *pict, int frame_index,
		int width, int height, fits *fit)
{
	/* when we pass a frame to the encoder, it may keep a reference to it
	 * internally; make sure we do not overwrite it here */
	if (av_frame_make_writable(pict) < 0)
		return 1;

	BYTE map[USHRT_MAX + 1];
	WORD tmp_pixel_value, hi, lo;
	int i;
	float pente;

	pente = computePente(&lo, &hi);
	
	for (i = 0; i <= USHRT_MAX; i++) {
		map[i] = round_to_BYTE((float) i * pente);
		if (map[i] == UCHAR_MAX)
			break;
	}
	if (i != USHRT_MAX + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= USHRT_MAX; i++)
			map[i] = UCHAR_MAX;
	}
	
	/* doing the WORD to BYTE conversion, bottom-up */
	if (fit->naxes[2] == 1) {
		int x, y;
		WORD *src = fit->pdata[RLAYER];
		BYTE *dst = pict->data[0];
		for (y = 0; y < fit->ry; y++) {
			int desty = fit->ry - y - 1;
			int srcpixel = y * fit->rx;
			int dstpixel = desty * fit->rx;
			for (x = 0; x < fit->rx; x++, srcpixel++, dstpixel++) {
				// linear scaling
				tmp_pixel_value = src[srcpixel] - lo;
				dst[dstpixel] = map[tmp_pixel_value];
			}
		}
	} else {
		int channel;
		for (channel = 0; channel < 3; channel++) {
			int x, y;
			WORD *src = fit->pdata[channel];
			BYTE *dst = pict->data[0];
			for (y = 0; y < fit->ry; y++) {
				int desty = fit->ry - y - 1;
				int srcpixel = y * fit->rx;
				int dstpixel = (desty * fit->rx * 3) + channel;
				for (x = 0; x < fit->rx; x++, srcpixel++, dstpixel+=3) {
					// linear scaling
					if (src[srcpixel] - lo < 0)
						tmp_pixel_value = 0;
					else 	tmp_pixel_value = src[srcpixel] - lo;
					dst[dstpixel] = map[tmp_pixel_value];
				}
			}
		}
	}
	return 0;
}

static AVFrame *get_video_frame(struct mp4_struct *ost, fits *input_image)
{
	AVCodecContext *c = ost->enc;

	/* check if we want to generate more frames -- why would we need that?
	 * BTW, STREAM_DURATION is invalid here, it should be long type
	 * if (av_compare_ts(ost->next_pts, c->time_base,
	 *         STREAM_DURATION, (AVRational){ 1, 1 }) >= 0)
	 *     return NULL;
	 */

	enum AVPixelFormat src_format =
		(input_image->naxes[2] == 1) ? AV_PIX_FMT_GRAY8 : AV_PIX_FMT_RGB24;

	/* if (target != input_image format) */
	if (c->pix_fmt != AV_PIX_FMT_RGB24) {
		if (!ost->sws_ctx) {
			ost->sws_ctx = sws_getContext(ost->src_w, ost->src_h, src_format,
					c->width, c->height, c->pix_fmt,
					SCALE_FLAGS, NULL, NULL, NULL);
			if (!ost->sws_ctx) {
				fprintf(stderr,
						"Could not initialize the conversion context\n");
				return NULL;
			}
		}
		fill_rgb_image(ost->tmp_frame, ost->next_pts, c->width, c->height, input_image);
		sws_scale(ost->sws_ctx,
				(const uint8_t * const *)ost->tmp_frame->data, ost->tmp_frame->linesize,
				0, ost->src_h, ost->frame->data, ost->frame->linesize);
	} else {
		fill_rgb_image(ost->frame, ost->next_pts, c->width, c->height, input_image);
	}

	ost->frame->pts = ost->next_pts++;

	return ost->frame;
}

/*
 * encode one video frame and send it to the muxer
 * return 1 when encoding is finished, 0 otherwise
 * https://ffmpeg.org/doxygen/3.1/group__lavc__encdec.html
 */
static int write_video_frame(struct mp4_struct *ost, fits *input_image)
{
	int ret;
	AVCodecContext *c = ost->enc;
	AVPacket pkt = { 0 };
	fprintf(stdout, "writing video frame\n");

	ost->frame = get_video_frame(ost, input_image);

	av_init_packet(&pkt);

	/* encode the image */
	ret = avcodec_send_frame(c, ost->frame);
	if (ret < 0) {
		fprintf(stderr, "Error encoding video frame: %s\n", av_err2str(ret));
		return 1;
	}

	ret = avcodec_receive_packet(c, &pkt);
	if (ret == 0) {
		ret = write_frame(ost->oc, &c->time_base, ost->st, &pkt);
		if (ret < 0) {
			fprintf(stderr, "Error while writing video frame: %s\n", av_err2str(ret));
			return 1;
		}
	}
	else if (ret == AVERROR(EINVAL)) {
		fprintf(stderr, "Error while getting video packet: %s\n", av_err2str(ret));
		return 1;
	}
	else if (ret == AVERROR_EOF) {
		fprintf(stderr, "End of stream met while adding a frame\n");
		return 1;
	}
	else if (ret == AVERROR(EAGAIN))
		return 0;
	//else retry later
	return 0;
}

static void flush_stream(struct mp4_struct *ost)
{
	int ret;
	AVPacket pkt = { 0 };
	av_init_packet(&pkt);

	ret = avcodec_send_frame(ost->enc, NULL);
	if (ret < 0) {
		fprintf(stderr, "Error encoding video frame: %s\n", av_err2str(ret));
		return;
	}

	do {
		ret = avcodec_receive_packet(ost->enc, &pkt);
		if (ret == 0) {
			ret = write_frame(ost->oc, &ost->enc->time_base, ost->st, &pkt);
			if (ret < 0) {
				fprintf(stderr, "Error while writing video frame: %s\n", av_err2str(ret));
				return;
			}
		}
		else if (ret == AVERROR(EINVAL)) {
			fprintf(stderr, "Error while getting video packet: %s\n", av_err2str(ret));
			return;
		}
	}
	while (ret != AVERROR_EOF) ;
}

static void close_stream(struct mp4_struct *ost)
{
	avcodec_free_context(&ost->enc);
	av_frame_free(&ost->frame);
	av_frame_free(&ost->tmp_frame);
	sws_freeContext(ost->sws_ctx);
	swr_free(&ost->swr_ctx);
}

/**************************************************************/
/* media file output */

struct mp4_struct *mp4_create(const char *filename, int dst_w, int dst_h, int fps, int nb_layers, int quality, int src_w, int src_h)
{
	struct mp4_struct *video_st;
	AVCodec *video_codec;
	int ret;
	AVDictionary *opt = NULL;

	if (filename == NULL || filename[0] == '\0' || dst_w%2 || dst_h%2 || fps <= 0 ||
			quality < 1 || quality > 5) {
		siril_log_message(_("Parameters for mp4 file creation were incorrect. Image dimension has to be a multiple of 2, fps and file name non nul, quality between 1 and 5."));
		return NULL;
	}

	/* Initialize libavcodec, and register all codecs and formats. */
	av_register_all();
	video_st = calloc(1, sizeof(struct mp4_struct));

	/* allocate the output media context */
	avformat_alloc_output_context2(&video_st->oc, NULL, NULL, filename);
	if (!video_st->oc) {
		fprintf(stdout, "Could not deduce output format from file extension: using mp4.\n");
		avformat_alloc_output_context2(&video_st->oc, NULL, "mp4", filename);
		if (!video_st->oc) {
			free(video_st);
			fprintf(stderr, "FFMPEG does not seem to support mp4 format, aborting.\n");
			return NULL;
		}
	}

	video_st->fmt = video_st->oc->oformat;
	/* disable unwanted features, is this the correct way? */
	video_st->fmt->audio_codec = AV_CODEC_ID_NONE;
	if (video_st->fmt->video_codec == AV_CODEC_ID_VP9) {
		/* force VP8 codec for webm, VP9 is not supported by Opera 12 */
		video_st->fmt->video_codec = AV_CODEC_ID_VP8;
	}
	video_st->src_w = src_w;
	video_st->src_h = src_h;
	video_st->bitrate = (quality + 1) * dst_w * dst_h / 2;

	/* Add the video stream and initialize the codecs. */
	if (video_st->fmt->video_codec != AV_CODEC_ID_NONE) {
		if (add_stream(video_st, &video_codec, video_st->fmt->video_codec, dst_w, dst_h, fps)) {
			avformat_free_context(video_st->oc);
			free(video_st);
			fprintf(stderr, "Could not add the video stream in the output film, aborting\n");
			return NULL;
		}
	} else {
		avformat_free_context(video_st->oc);
		fprintf(stderr, "Error setting video codec for output file.");
		return NULL;
	}

	/* Now that all the parameters are set, we can open the video codec
	 * and allocate the necessary encode buffers. */
	if (open_video(video_codec, video_st, opt, nb_layers)) {
		avformat_free_context(video_st->oc);
		avcodec_free_context(&video_st->enc);
		free(video_st);
		return NULL;
	}

	av_dump_format(video_st->oc, 0, filename, 1);

	/* open the output file, if needed */
	if (!(video_st->fmt->flags & AVFMT_NOFILE)) {
		if ((ret = avio_open(&video_st->oc->pb, filename, AVIO_FLAG_WRITE)) < 0) {
			fprintf(stderr, "Could not open '%s': %s\n", filename, av_err2str(ret));
			avformat_free_context(video_st->oc);
			avcodec_free_context(&video_st->enc);
			av_frame_free(&video_st->frame);
			if (video_st->tmp_frame) av_frame_free(&video_st->tmp_frame);
			free(video_st);
			return NULL;
		}
	}

	/* Write the stream header, if any. */
	if ((ret = avformat_write_header(video_st->oc, &opt)) < 0) {
		fprintf(stderr, "Error occurred when opening output file: %s\n", av_err2str(ret));
		avformat_free_context(video_st->oc);
		avcodec_free_context(&video_st->enc);
		av_frame_free(&video_st->frame);
		if (video_st->tmp_frame) av_frame_free(&video_st->tmp_frame);
		free(video_st);
		return NULL;
	}

	return video_st;
}

int mp4_add_frame(struct mp4_struct *video_st, fits *image) {
	/* select the stream to encode */
	if (av_compare_ts(video_st->next_pts, video_st->enc->time_base,
				video_st->next_pts, video_st->enc->time_base) <= 0) {
		return write_video_frame(video_st, image);
	}

	/* should never happen, because we don't limit it */
	fprintf(stderr, "End of video stream\n");
	return -1;
}

int mp4_close(struct mp4_struct *video_st) {
	flush_stream(video_st);

	/* Write the trailer, if any. The trailer must be written before you
	 * close the CodecContexts open when you wrote the header; otherwise
	 * av_write_trailer() may try to use memory that was freed on
	 * av_codec_close(). */
	av_write_trailer(video_st->oc);

	/* Close each codec. */
	close_stream(video_st);

	if (!(video_st->fmt->flags & AVFMT_NOFILE))
		/* Close the output file. */
		avio_closep(&video_st->oc->pb);

	/* free the stream */
	avformat_free_context(video_st->oc);
	avcodec_free_context(&video_st->enc);
	av_frame_free(&video_st->frame);
	if (video_st->tmp_frame)
		av_frame_free(&video_st->tmp_frame);

	return 0;
}

#endif
