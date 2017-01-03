/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gui/callbacks.h"
#include "io/films.h"
#include "core/proto.h"		// for fits_flip_top_to_bottom()

static int pixfmt_gray, pixfmt_rgb, pixfmt_gray16, pixfmt_rgb48;

supported_film_list supported_film[] = {
	{"avi"},
	{"mpg"},
	{"mpeg"},
	{"mov"},
	{"mp4"},
	{"webm"}
};

int get_nb_film_ext_supported() {
	return sizeof(supported_film) / sizeof(supported_film_list);
}

static void film_init_struct(struct film_struct *film) {
	memset(film, 0, sizeof(struct film_struct));
}

/* Check different film extensions supported in supported_film[].
 * The function return 0 if an extension is found, 1 else */
int check_for_film_extensions(const char *extension) {
	int i, nb_film;

	nb_film = get_nb_film_ext_supported();
	for (i = 0; i < nb_film; i++) {
		if (!g_ascii_strcasecmp(extension, supported_film[i].extension)) return 0;
	}
	return 1;
}

int film_open_file(const char *sourcefile, struct film_struct *film) {
	film_init_struct(film);
	/* Initialize the library itself. */
	FFMS_Init(0, 0);

	/* Index the source file. Note that this example does not index any audio tracks. */
	film->errmsg = malloc(FILM_ERROR_LENGTH);
	film->errinfo.Buffer      = film->errmsg;
	film->errinfo.BufferSize  = FILM_ERROR_LENGTH;
	film->errinfo.ErrorType   = FFMS_ERROR_SUCCESS;
	film->errinfo.SubType     = FFMS_ERROR_SUCCESS;

	/* First, try to read it, in case the film has already been indexed and the index is
	 * present on disk.
	 * Index is saved using this pattern for the file name "film_filename.idx" */
	FFMS_Index *index;
#ifdef HAVE_FFMS2_2
	FFMS_Indexer *indexer;
#endif
	char *idxfilename;
	idxfilename = malloc(strlen(sourcefile) + 5);
	sprintf(idxfilename, "%s.idx", sourcefile);
	index = FFMS_ReadIndex(idxfilename, &film->errinfo);
	if (index == NULL) {
#ifdef HAVE_FFMS2_2
		/* we need to create the indexer */
		indexer = FFMS_CreateIndexer(sourcefile, &film->errinfo);
		if (indexer == NULL) {
#else
		/* we need to create the index */
		index = FFMS_MakeIndex(sourcefile, 0, 0, NULL, NULL,
				FFMS_IEH_ABORT, NULL, NULL, &film->errinfo);
		if (index == NULL) {
#endif
			/* handle error (print errinfo.Buffer somewhere) */
			fprintf(stderr, "FILM error: %s\n", film->errmsg);
			free(idxfilename);
			return FILM_ERROR;
		}
#ifdef HAVE_FFMS2_2
		/* we need to create the index */
		index = FFMS_DoIndexing2(indexer, FFMS_IEH_ABORT, &film->errinfo);
		if (index == NULL) {
			/* handle error (print errinfo.Buffer somewhere) */
			fprintf(stderr, "FILM error: %s\n", film->errmsg);
			free(idxfilename);
			return FILM_ERROR;
		}
#endif

		/* write the index for future openings */
		if (FFMS_WriteIndex(idxfilename, index, &film->errinfo)) {
			fprintf(stderr, "FILM: Could not save index file: %s\n", film->errmsg);
		}
		else fprintf(stdout, "FILM: index saved into file '%s'\n", idxfilename);
	} else fprintf(stdout, "FILM: loaded previously computed index from file '%s'\n", idxfilename);
	free(idxfilename);

	/* Retrieve the track number of the first video track */
	int trackno = FFMS_GetFirstTrackOfType(index, FFMS_TYPE_VIDEO, &film->errinfo);
	if (trackno < 0) {
		/* no video tracks found in the file, this is bad and you should handle it */
		/* (print the errmsg somewhere) */
		fprintf(stderr, "FILM error: %s\n", film->errmsg);
		return FILM_ERROR;
	}

	/* We now have enough information to create the video source object */
	film->videosource = FFMS_CreateVideoSource(sourcefile, trackno, index, 1, FFMS_SEEK_NORMAL, &film->errinfo);
	if (film->videosource == NULL) {
		/* handle error (you should know what to do by now) */
		fprintf(stderr, "FILM error: %s\n", film->errmsg);
		return FILM_ERROR;
	}

	/* Since the index is copied into the video source object upon its creation,
	we can and should now destroy the index object. */
	FFMS_DestroyIndex(index);

	/* Retrieve video properties so we know what we're getting.
	As the lack of the errmsg parameter indicates, this function cannot fail. */
	const FFMS_VideoProperties *videoprops = FFMS_GetVideoProperties(film->videosource);

	/* Now you may want to do something with the info, like check how many frames the video has */
	film->frame_count = videoprops->NumFrames;

	/* Get the first frame for examination so we know what we're getting. This is required
	because resolution and colorspace is a per frame property and NOT global for the video. */
	const FFMS_Frame *propframe = FFMS_GetFrame(film->videosource, 0, &film->errinfo);

	/* Now you may want to do something with the info; particularly interesting values are:
	propframe->EncodedWidth; (frame width in pixels)
	propframe->EncodedHeight; (frame height in pixels)
	propframe->EncodedPixelFormat; (actual frame colorspace)
	*/
	film->width = propframe->EncodedWidth;
	film->height = propframe->EncodedHeight;
	/* pixel format, giving the number of layers, is guessed here from the original format.
	 * However, a film containing gray images can be encoded as RGB32, in order to keep
	 * the same format between the black and white and the color version of a camera. This
	 * has to be detected when reading the frames, checking values in the layers, lame and
	 * costly. */
	pixfmt_gray = FFMS_GetPixFmt("gray8");		// return value is PIX_FMT_NONE
	pixfmt_rgb = FFMS_GetPixFmt("rgb24");
	pixfmt_gray16 = FFMS_GetPixFmt("gray16");
	pixfmt_rgb48 = FFMS_GetPixFmt("rgb48");
	if (propframe->EncodedPixelFormat == pixfmt_gray16 ||
			propframe->EncodedPixelFormat == pixfmt_rgb48) {
		film->nb_layers = 0;
		film->pixfmt = 0;
		fprintf(stderr, "FILM: 16-bit pixel depth films are not supported yet.\n");
		return FILM_ERROR;
	}
	else if (propframe->EncodedPixelFormat == pixfmt_gray) {
		film->nb_layers = 1;
		film->pixfmt = pixfmt_gray;
	} else {
		film->nb_layers = 3;
		film->pixfmt = pixfmt_rgb;
	}

	/* If you want to change the output colorspace or resize the output frame size,
	now is the time to do it. IMPORTANT: This step is also required to prevent
	resolution and colorspace changes midstream. You can always tell a frame's
	original properties by examining the Encoded properties in FFMS_Frame. */
	/* See libavutil/pixfmt.h for the list of pixel formats/colorspaces.
	To get the name of a given pixel format, strip the leading PIX_FMT_
	and convert to lowercase. For example, PIX_FMT_YUV420P becomes "yuv420p". */

	/* A -1 terminated list of the acceptable output formats. */
	int pixfmts[2];
	pixfmts[0] = film->pixfmt;
	pixfmts[1] = -1;

	if (FFMS_SetOutputFormatV2(film->videosource, pixfmts,
			propframe->EncodedWidth, propframe->EncodedHeight,
			FFMS_RESIZER_BICUBIC, &film->errinfo)) {
		/* handle error */
		fprintf(stderr, "FILM error: %s\n", film->errmsg);
		return FILM_ERROR;
	}

	film->filename = strdup(sourcefile);
	fprintf(stdout, "FILM: successfully opened the video file %s, %d frames\n",
			film->filename, film->frame_count);
	return FILM_SUCCESS;
}

static int randPixel(int nb_pixels) {
	return rand() % nb_pixels;
}

static int *randomIndex(int n) {
	srand(time(NULL));
	int *index = malloc(n * sizeof (int));
	int i = 0;
	int x = 0;
	int tmp = 0;

	/* make index */
	for (i = 0; i < n; i++) {
		index[i] = i;
	}

	/* mix index */
	for (i = 0; i < n; i++) {
		x = randPixel(n);
		tmp = index[i];
		index[i] = index[x];
		index[x] = tmp;
	}	return index;
}

int film_read_frame(struct film_struct *film, int frame_no, fits *fit) {
	/* now we're ready to actually retrieve the video frames */
	int nb_pixels, convert_rgb_to_gray = 0;

	if (film->videosource == 0x00) {
		siril_log_message(_("FILM ERROR: incompatible format\n"));
		return FILM_ERROR;
	}
	const FFMS_Frame *frame = FFMS_GetFrame(film->videosource, frame_no, &film->errinfo);
	WORD *ptr;
	if (frame == NULL) {
		/* handle error */
		fprintf(stderr, "FILM error: %s\n", film->errmsg);
		return FILM_ERROR;
	}

	nb_pixels = film->width * film->height;

	/* detect gray images encoded in RGB24 movies */
	if (frame->ConvertedPixelFormat == pixfmt_rgb) {
		// browse random pixels in the first frame (100 pixels) and check
		// value in two or three layers
		// convert_rgb_to_gray = 1;
		if (frame_no == 0) {
			int pixel_tested = 0, n = 0;
			int *randIndex = randomIndex(nb_pixels * 3);
			do {
				int px = randIndex[n];
				px = px - (px % 3);
				++n;
				WORD r, g, b;

				r = frame->Data[0][px + 0];
				g = frame->Data[0][px + 1];
				b = frame->Data[0][px + 2];

				if (r == g && r == b && g == b) {
					convert_rgb_to_gray = 1;
				} else {
					/* no gray image, we can go out of the loop */
					convert_rgb_to_gray = 0;
					break;
				}
				// we reject pure black and pure white in the comparison
				if (r != 0 && r != 255)
					++pixel_tested;

			} while (pixel_tested < 100 && n < nb_pixels * 3);
			free(randIndex);
			printf("total n = %d et k = %d et npixel = %d\n", n, pixel_tested, nb_pixels);
		}
	}
	if (convert_rgb_to_gray)
		film->nb_layers = 1;

	/* do something with frame */

	if ((ptr = realloc(fit->data, nb_pixels * film->nb_layers * sizeof(WORD)))
			== NULL) {
		fprintf(stderr,"FILM: NULL realloc for FITS data\n");
		free(fit->data);
		return -1;
	}
	memset(fit, 0, sizeof(fits));
	fit->data = ptr;
	fit->naxes[0] = fit->rx = film->width;
	fit->naxes[1] = fit->ry = film->height;
	fit->naxes[2] = film->nb_layers;
	if (fit->naxes[2] == 1) {
		fit->naxis = 2;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data;
		fit->pdata[BLAYER] = fit->data;
	} else {
		fit->naxis = 3;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + nb_pixels;
		fit->pdata[BLAYER] = fit->data + nb_pixels * 2;
	}
	fit->bitpix = BYTE_IMG;
	//fit->mini = 0;
	//fit->maxi = 255;
	/* putting this above also requires the max[*] to be = 255. Besides, this overrides the
	 * default min/max behavior of Siril. */

	if (frame->ConvertedPixelFormat == pixfmt_gray || convert_rgb_to_gray) {
		int i, j;
		for (i = 0, j = 0; i < nb_pixels; i++) {
			fit->pdata[RLAYER][i] = frame->Data[0][j];
			if (convert_rgb_to_gray)
				j += 3;
			else
				++j;
		}
	} else if (frame->ConvertedPixelFormat == pixfmt_rgb) {
		int i, j;
		for (i = 0, j = 0; i < nb_pixels; i++) {
			fit->pdata[RLAYER][i] = frame->Data[0][j++];
			fit->pdata[GLAYER][i] = frame->Data[0][j++];
			fit->pdata[BLAYER][i] = frame->Data[0][j++];
		}
	}
	else {
		// format is not one we set, should happen only if return value of the file
		// opening was not used to discard the file
		fprintf(stderr, "FILM: format not understood\n");
		return FILM_ERROR;
	}
	fits_flip_top_to_bottom(fit);

	return FILM_SUCCESS;
}

void film_close_file(struct film_struct *film) {
	/* now it's time to clean up */
	FFMS_DestroyVideoSource(film->videosource);
}

void film_display_info(struct film_struct *film) {
	fprintf(stdout, "\n============= FILM file info =============\n");
	fprintf(stdout, "file name: %s\n", film->filename);
//	fprintf(stdout, "video stream number: %d\n", film->stream_number);
	fprintf(stdout, "image size: %d x %d\n", film->width, film->height);
	fprintf(stdout, "number of layers: %d\n", film->nb_layers);
	fprintf(stdout, "frame count: %d\n", film->frame_count);
	fprintf(stdout, "==========================================\n");
}

#endif
