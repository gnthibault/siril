#ifndef _FILMS_H_
#define _FILMS_H_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)

#include <ffms.h>

#define FILM_SUCCESS 0
#define FILM_ERROR -1

#define FILM_ERROR_LENGTH 300

typedef struct {
	char *extension;
}supported_film_list;

extern supported_film_list supported_film[];	//supported film extensions

struct film_struct {
	FFMS_VideoSource *videosource;
	FFMS_ErrorInfo errinfo;
	int pixfmt;
	char *errmsg;

	int width, height;
	int nb_layers;		// 1 for gray, 3 for rgb, 0 for uninit
	int frame_count;

	char *filename;
};

/* external functions */
int get_nb_film_ext_supported();
int check_for_film_extensions(const char *extension);
int film_open_file(const char *sourcefile, struct film_struct *film);
void film_close_file(struct film_struct *film);
int film_read_frame(struct film_struct *film, int frame_no, fits *fit);
void film_display_info(struct film_struct *film);

#endif
#endif

