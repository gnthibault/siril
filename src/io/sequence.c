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

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <dirent.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <libgen.h>

#include "core/siril.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/plot.h"
#include "io/ser.h"
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
#include "io/films.h"
#endif
#ifdef HAVE_OPENCV
#include "opencv/opencv.h"
#endif
#ifdef HAVE_FFMPEG
#include "io/mp4_output.h"
#endif
#include "io/avi_pipp/avi_writer.h"
#include "io/single_image.h"
#include "gui/histogram.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"	// clear_stars_list
#include "algos/quality.h"
#include "registration/registration.h"	// for update_reg_interface
#include "stacking/stacking.h"	// for update_stack_interface

static void fillSeqAviExport() {
	char width[6], height[6], fps[7];
	GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));
	GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));

	g_snprintf(width, sizeof(width), "%d", com.seq.rx);
	g_snprintf(height, sizeof(width), "%d", com.seq.ry);
	gtk_entry_set_text(widthEntry, width);
	gtk_entry_set_text(heightEntry, height);
	if (com.seq.type == SEQ_SER) {
		GtkEntry *entryAviFps = GTK_ENTRY(lookup_widget("entryAviFps"));

		if (com.seq.ser_file != NULL) {
			if (com.seq.ser_file->fps <= 0.0) {
				g_snprintf(fps, sizeof(fps), "25.000");
			} else {
				g_snprintf(fps, sizeof(fps), "%2.3lf", com.seq.ser_file->fps);
			}
			gtk_entry_set_text(entryAviFps, fps);
		}
	}
}

/* when opening a file outside the main sequence loading system and that file
 * is a sequence (SER/AVI), this function is called to load this sequence. */
int read_single_sequence(char *realname, int imagetype) {
	int retval=3;		// needs to return 3 if ok !!!
	char *name = strdup(realname);
	gchar *dirname = g_path_get_dirname(realname);
	if (!changedir(dirname))
		writeinitfile();
	g_free(dirname);

	if (check_only_one_film_seq(realname)) retval = 1;
	else {
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	const char *ext;
#endif
		switch (imagetype) {
			case TYPESER:
				name[strlen(name)-1] = 'q';
				break;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
			case TYPEAVI:
				ext = get_filename_ext(realname);
				assert(ext);
				int len = strlen(ext);
				strncpy(name+strlen(name)-len, "seq", len);
				break;
#endif
			default:
				retval = 1;
		}
		gchar *fname = g_path_get_basename(name);
		if (!set_seq(fname)) {
			/* if it loads, make it selected and only element in the list of sequences */
			control_window_switch_to_tab(IMAGE_SEQ);
			GtkComboBoxText *combo_box_text = GTK_COMBO_BOX_TEXT(lookup_widget("sequence_list_combobox"));
			gtk_combo_box_text_remove_all(combo_box_text);
			gchar *rname = g_path_get_basename(realname);
			gtk_combo_box_text_append(combo_box_text, 0, rname);
			g_signal_handlers_block_by_func(GTK_COMBO_BOX(combo_box_text), on_seqproc_entry_changed, NULL);
			gtk_combo_box_set_active(GTK_COMBO_BOX(combo_box_text), 0);
			g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(combo_box_text), on_seqproc_entry_changed, NULL);
			g_free(rname);
		}
		else retval = 1;
		g_free(fname);
	}
	free(name);
	return retval;
}

/* Find sequences in CWD and create .seq files.
 * In the current working directory, looks for sequences of fits files or files
 * already representing sequences like SER and AVI formats and builds the
 * corresponding sequence files.
 * Called when changing wd with name == NULL or when an explicit root name is
 * given in the GUI or when searching for sequences.
 */
int check_seq(int force) {
	char *basename;
	int curidx, fixed;
	DIR *dir;
	struct dirent *file;
	sequence **sequences;
	int i, nb_seq = 0, max_seq = 10;

	if (!com.wd) {
		siril_log_message(_("Current working directory is not set, aborting.\n"));
		return 1;
	}
	if ((dir = opendir(com.wd)) == NULL) {
		fprintf(stderr, "working directory cannot be opened.\n");
		free(com.wd);
		com.wd = NULL;
		return 1;
	}

	sequences = malloc(sizeof(sequence *) * max_seq);

	while ((file = readdir(dir)) != NULL) {
		sequence *new_seq;
		int fnlen = strlen(file->d_name);
		if (fnlen < 4) continue;
		const char *ext = get_filename_ext(file->d_name);
		if (!ext) continue;
		if (!strcasecmp(ext, "ser")) {
			struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
			ser_init_struct(ser_file);
			if (ser_open_file(file->d_name, ser_file))
				continue;
			new_seq = calloc(1, sizeof(sequence));
			initialize_sequence(new_seq, TRUE);
			new_seq->seqname = g_strndup(file->d_name, fnlen-4);
			new_seq->beg = 0;
			new_seq->end = ser_file->frame_count-1;
			new_seq->number = ser_file->frame_count;
			new_seq->type = SEQ_SER;
			new_seq->ser_file = ser_file;
			sequences[nb_seq] = new_seq;
			nb_seq++;
			fprintf(stdout, "Found a SER sequence (number %d)\n", nb_seq);
			set_progress_bar_data(NULL, PROGRESS_PULSATE);
		}
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		else if (!check_for_film_extensions(ext)) {
			struct film_struct *film_file = malloc(sizeof(struct film_struct));
			if (film_open_file(file->d_name, film_file)) {
				free(film_file);
				continue;
			}
			new_seq = calloc(1, sizeof(sequence));
			initialize_sequence(new_seq, TRUE);
			int len = strlen(ext);
			new_seq->seqname = g_strndup(file->d_name, fnlen-(len+1));
			new_seq->beg = 0;
			new_seq->end = film_file->frame_count-1;
			new_seq->number = film_file->frame_count;
			new_seq->type = SEQ_AVI;
			new_seq->film_file = film_file;
			sequences[nb_seq] = new_seq;
			nb_seq++;
			fprintf(stdout, "Found a AVI sequence (number %d)\n", nb_seq);
			set_progress_bar_data(NULL, PROGRESS_PULSATE);
		}
#endif

		else if (!strcasecmp(ext, com.ext+1)) {
			if (!get_index_and_basename(file->d_name, &basename, &curidx, &fixed)) {
				int current_seq = -1;
				/* search in known sequences if we already have it */
				for (i=0; i<nb_seq; i++) {
					if (!strcmp(sequences[i]->seqname, basename)) {
						current_seq = i;
					}
				}
				/* not found */
				if (current_seq == -1) {
					new_seq = calloc(1, sizeof(sequence));
					initialize_sequence(new_seq, TRUE);
					new_seq->seqname = basename;
					new_seq->beg = INT_MAX;
					new_seq->end = 0;
					new_seq->fixed = fixed;
					sequences[nb_seq] = new_seq;
					current_seq = nb_seq;
					nb_seq++;
					fprintf(stdout, "Found a sequence (number %d) with base name"
							" \"%s\", looking for first and last indexes.\n",
							nb_seq, basename);
					set_progress_bar_data(NULL, PROGRESS_PULSATE);
				}
				if (curidx < sequences[current_seq]->beg)
					sequences[current_seq]->beg = curidx;
				if (curidx > sequences[current_seq]->end)
					sequences[current_seq]->end = curidx;
				if (fixed > sequences[current_seq]->fixed)
					sequences[current_seq]->fixed = fixed;
			}
		}
		if (nb_seq == max_seq) {
			max_seq *= 2;
			sequence **tmp = realloc(sequences, sizeof(sequence *) * max_seq);
			if (tmp)
				sequences = tmp;
			else {
				siril_log_message(_("Could not allocate more space for the large number of sequences found.\n"));
				break;
			}
		}
	}
	closedir(dir);
	if (nb_seq > 0) {
		int retval = 1;
		for (i=0; i<nb_seq; i++) {
			if (sequences[i]->beg != sequences[i]->end) {
				char msg[200];
				sprintf(msg, _("sequence %d, found: %d to %d"),
						i+1, sequences[i]->beg, sequences[i]->end);
				set_progress_bar_data(msg, PROGRESS_NONE);
				if (!buildseqfile(sequences[i], force) && retval)
					retval = 0;	// at least one succeeded to be created
			}
			free_sequence(sequences[i], TRUE);
		}
		free(sequences);
		return retval;
	}
	free(sequences);
	return 1;	// no sequence found
}

/* Check for on film sequence of the name passed in arguement
 * Returns 0 if OK */
int check_only_one_film_seq(char* name) {
	int retval = 1;
	DIR *dir;
	sequence *new_seq = NULL;

	if (!com.wd) {
		siril_log_message(_("Current working directory is not set, aborting.\n"));
		return 1;
	}
	if ((dir = opendir(com.wd)) == NULL) {
		fprintf(stderr, "working directory cannot be opened.\n");
		free(com.wd);
		com.wd = NULL;
		return 1;
	}

	int fnlen = strlen(name);
	const char *ext = get_filename_ext(name);

	if (!strcasecmp(ext, "ser")) {
		struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
		ser_init_struct(ser_file);
		if (ser_open_file(name, ser_file)) {
			closedir(dir);
			return 1;
		}

		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		new_seq->seqname = g_strndup(name, fnlen-4);
		new_seq->beg = 0;
		new_seq->end = ser_file->frame_count-1;
		new_seq->number = ser_file->frame_count;
		new_seq->type = SEQ_SER;
		new_seq->ser_file = ser_file;
	}
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	else if (!check_for_film_extensions(ext)) {
		struct film_struct *film_file = malloc(sizeof(struct film_struct));
		if (film_open_file(name, film_file)) {
			free(film_file);
			closedir(dir);
			return 1;
		}
		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		int len = strlen(ext);
		new_seq->seqname = g_strndup(name, fnlen-len-1);
		new_seq->beg = 0;
		new_seq->end = film_file->frame_count-1;
		new_seq->number = film_file->frame_count;
		new_seq->type = SEQ_AVI;
		new_seq->film_file = film_file;
		fprintf(stdout, "Found a AVI sequence\n");
	}
#endif
	closedir(dir);
	if (!new_seq) return 0;
	if (new_seq->beg != new_seq->end) {
		if (!buildseqfile(new_seq, 0) && retval)
			retval = 0;
	}
	free_sequence(new_seq, TRUE);
	return retval;
}

/* load a sequence and initializes everything that relates */
int set_seq(const char *name){
	sequence *seq;
	int image_to_load;
	char *basename;
	
	if ((seq = readseqfile(name)) == NULL) {
		fprintf(stderr, "could not load sequence %s\n", name);
		return 1;
	}
	free_image_data();

	image_to_load = sequence_find_refimage(seq);
	if (seq_read_frame(seq, image_to_load, &gfit)) {
		fprintf(stderr, "could not load first image from sequence\n");
		free(seq);
		return 1;
	}

	/* initialize sequence-related runtime data */
	seq->rx = gfit.rx; seq->ry = gfit.ry;
	seq->current = image_to_load;

	if (seq->nb_layers == -1 || seq->nb_layers != gfit.naxes[2]) {	// not init yet, first loading of the sequence
		seq->nb_layers = gfit.naxes[2];
		seq->regparam = calloc(seq->nb_layers, sizeof(regdata *));
		seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
		writeseqfile(seq);
	}

	basename = g_path_get_basename(seq->seqname);
	siril_log_message(_("Sequence loaded: %s (%d->%d)\n"), basename, seq->beg,
			seq->end);
	g_free(basename);
	/* Sequence is stored in com.seq for now */
	free_sequence(&com.seq, FALSE);
	memcpy(&com.seq, seq, sizeof(sequence));

	if (seq->nb_layers > 1)
		show_rgb_window();
	else hide_rgb_window();
	init_layers_hi_and_lo_values(MIPSLOHI); // set some hi and lo values in seq->layers,
	set_cutoff_sliders_max_values();// update min and max values for contrast sliders
	set_cutoff_sliders_values();	// update values for contrast sliders for this image
	seqsetnum(image_to_load);	// set limits for spin button and display loaded filenum
	set_layers_for_assign();	// set default layers assign and populate combo box
	set_layers_for_registration();	// set layers in the combo box for registration
	fill_sequence_list(seq, 0);	// display list of files in the sequence
	set_output_filename_to_sequence_name();
	sliders_mode_set_state(com.sliders);
	initialize_display_mode();

	/* initialize image-related runtime data */
	set_display_mode();		// display the display mode in the combo box
	display_filename();		// display filename in gray window
	adjust_exclude(image_to_load, FALSE);	// check or uncheck excluded checkbox
	adjust_refimage(image_to_load);	// check or uncheck reference image checkbox
	set_prepro_button_sensitiveness(); // enable or not the preprobutton
	update_reg_interface(FALSE);	// change the registration prereq message
	update_stack_interface(FALSE);	// get stacking info and enable the Go button
	adjust_reginfo();		// change registration displayed/editable values
	update_gfit_histogram_if_needed();
	adjust_sellabel();
	fillSeqAviExport();	// fill GtkEntry of export box

	/* update menus */
	update_MenuItem();
	/* update parameters */
	set_GUI_CAMERA();
	set_GUI_photometry();

	/* redraw and display image */
	show_main_gray_window();
	close_tab();	//close Green and Blue Tab if a 1-layer sequence is loaded
	adjust_vport_size_to_image();	// resize viewports to the displayed image size
	redraw(com.cvport, REMAP_ALL);
	drawPlot();

	update_used_memory();
	return 0;
}

/* Load image number index from the sequence and display it.
 * if load_it is true, dest is assumed to be gfit
 * TODO: cut that method in two, with an internal func taking a filename and a fits
 */
int seq_load_image(sequence *seq, int index, fits *dest, gboolean load_it) {
	seq->current = index;
	clear_stars_list();
	clear_histograms();
	gfit.maxi = 0;
	// what else needs to be cleaned?
	if (load_it) {
		set_cursor_waiting(TRUE);
		if (seq_read_frame(seq, index, dest)) {
			set_cursor_waiting(FALSE);
			return 1;
		}
		set_fwhm_star_as_star_list(seq);// display the fwhm star if possible
		if (com.sliders != USER) {
			init_layers_hi_and_lo_values(com.sliders);
			sliders_mode_set_state(com.sliders);
			set_cutoff_sliders_max_values();// update min and max values for contrast sliders
			set_cutoff_sliders_values();	// update values for contrast sliders for this image
			set_display_mode();		// display the display mode in the combo box
		}
		if (copy_rendering_settings_when_chained(TRUE))
			redraw(com.cvport, REMAP_ALL);
		else
			redraw(com.cvport, REMAP_ONLY);
		redraw_previews();		// redraw registration preview areas
		display_filename();		// display filename in gray window
		adjust_reginfo();		// change registration displayed/editable values
		calculate_fwhm(com.vport[com.cvport]);
		update_gfit_histogram_if_needed();
		set_cursor_waiting(FALSE);
	}
	/* change the displayed value in the spin button to have the real file number
	 * instead of the index of the adjustment */

	undo_flush();
	/* initialize menu gui */
	update_MenuItem();

	display_image_number(index);
	sequence_list_change_current();
	adjust_exclude(index, FALSE);	// check or uncheck excluded checkbox
	adjust_refimage(index);	// check or uncheck reference image checkbox
	update_used_memory();
	return 0;
}

/*****************************************************************************
 *              SEQUENCE FUNCTIONS FOR NON-OPENED SEQUENCES                  *
 * **************************************************************************/

/* Get the filename of an image in a sequence.
 * Return value is the same as the name_buf argument, which must be
 * pre-allocated to at least 256 characters. If sequence has no file names, a
 * description like image "42 from awesome_mars.ser" is made. */
char *seq_get_image_filename(sequence *seq, int index, char *name_buf) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return fit_sequence_get_image_filename(seq, index, name_buf, TRUE);
		case SEQ_SER:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			snprintf(name_buf, 255, _("%d from %s.ser"), index, seq->seqname);
			name_buf[255] = '\0';
			return name_buf;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		case SEQ_AVI:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			snprintf(name_buf, 255, "%s_%d", seq->seqname, index);
			name_buf[255] = '\0';
			return name_buf;
#endif
		case SEQ_INTERNAL:
			snprintf(name_buf, 255, "%s_%d", seq->seqname, index);
			name_buf[255] = '\0';
			return name_buf;
	}
	return NULL;
}

/* Read an entire image from a sequence, inside a pre-allocated fits.
 * Opens the file, reads data, closes the file.
 */
int seq_read_frame(sequence *seq, int index, fits *dest) {
	char filename[256];
	assert(index < seq->number);
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			if (readfits(filename, dest, NULL)) {
				siril_log_message(_("Could not load image %d from sequence %s\n"),
						index, seq->seqname); 
				return 1;
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			if (ser_read_frame(seq->ser_file, index, dest)) {
				siril_log_message(_("Could not load frame %d from SER sequence %s\n"),
						index, seq->seqname); 
				return 1;
			}
			break;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		case SEQ_AVI:
			assert(seq->film_file);
			if (film_read_frame(seq->film_file, index, dest)) {
				siril_log_message(_("Could not load frame %d from AVI sequence %s\n"),
						index, seq->seqname); 
				return 1;
			}
			// should dest->maxi be set to 255 here?
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			copyfits(seq->internal_fits[index], dest, CP_FORMAT, -1);
			dest->data = seq->internal_fits[index]->data;
			dest->pdata[0] = seq->internal_fits[index]->pdata[0];
			dest->pdata[1] = seq->internal_fits[index]->pdata[1];
			dest->pdata[2] = seq->internal_fits[index]->pdata[2];
			break;
	}
	image_find_minmax(dest, 0);
	return 0;
}

/* same as seq_read_frame above, but creates an image the size of the selection
 * rectangle only. layer is set to the layer number in the read partial frame.
 * The partial image result is only one-channel deep, so it cannot be used to
 * have a partial RGB image. */
int seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry) {
	char filename[256];
	fits tmp_fit;
	memset(&tmp_fit, 0, sizeof(fits));
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			if (readfits_partial(filename, layer, dest, area, do_photometry)) {
				siril_log_message(_("Could not load partial image %d from sequence %s\n"),
						index, seq->seqname); 
				return 1;
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			/* TODO: build a FITS from ser_read_opened_partial() */
			if (ser_read_frame(seq->ser_file, index, &tmp_fit)) {
				siril_log_message(_("Could not load frame %d from SER sequence %s\n"),
						index, seq->seqname); 
				return 1;
			}
			extract_region_from_fits(&tmp_fit, layer, dest, area);
			clearfits(&tmp_fit);
			break;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		case SEQ_AVI:
			assert(seq->film_file);
			if (film_read_frame(seq->film_file, index, &tmp_fit)) {
				siril_log_message(_("Could not load frame %d from AVI sequence %s\n"),
						index, seq->seqname); 
				return 1;
			}
			extract_region_from_fits(&tmp_fit, layer, dest, area);
			clearfits(&tmp_fit);
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			extract_region_from_fits(seq->internal_fits[index], 0, dest, area);
			break;
	}
	return 0;
}

/*****************************************************************************
 *                 SEQUENCE FUNCTIONS FOR OPENED SEQUENCES                   *
 * **************************************************************************/

/* locks cannot be probed to see if they are init or not, so we have to keep
 * all of them in the same state, which is initialized if the array is non-nul. */
static int _allocate_sequence_locks(sequence *seq) {
#ifdef _OPENMP
	if (!seq->fd_lock) {
		int i;
		seq->fd_lock = malloc(seq->number * sizeof(omp_lock_t));
		if (!seq->fd_lock) {
			fprintf(stderr, "Allocation error when opening images, aborting\n");
			return 1;
		}

		for (i=0; i<seq->number; i++)
			omp_init_lock(&seq->fd_lock[i]);
	}
#endif
	return 0;
}

/* open image for future intensive operations (read only) */
int seq_open_image(sequence *seq, int index) {
	int status = 0;
	char filename[256];
	switch (seq->type) {
		case SEQ_REGULAR:
			if (!seq->fptr) {
				seq->fptr = calloc(seq->number, sizeof(fitsfile *));
				if (!seq->fptr) {
				       fprintf(stderr, "Allocation error when opening images, aborting\n");
			       	       return 1;
				}
			}
			if (_allocate_sequence_locks(seq))
				return 1;

			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			fits_open_diskfile(&seq->fptr[index], filename, READONLY, &status);
			if (status) {
				fits_report_error(stderr, status);
				return status;
			}
			/* should we check image parameters here? such as bitpix or naxis */
			break;
		case SEQ_SER:
			assert(seq->ser_file->fd > 0);
			break;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		case SEQ_AVI:
			siril_log_message(_("This operation is not supported on AVI sequences (seq_open_image)\n"));
			return 1;
#endif
		case SEQ_INTERNAL:
			siril_log_message(_("This operation is not supported on internal sequences (seq_open_image)\n"));
			return 1;
	}
	return 0;
}

/* close opened images, only useful for regular FITS sequences */
void seq_close_image(sequence *seq, int index) {
	int status = 0;
	switch (seq->type) {
		case SEQ_REGULAR:
			if (seq->fptr && seq->fptr[index]) {
				fits_close_file(seq->fptr[index], &status);
				seq->fptr[index] = NULL;
			}
			break;
		default:
			break;
	}
}

/* read a region in a layer of an opened file from a sequence.
 * The buffer must have been allocated to the size of the area. */
int seq_opened_read_region(sequence *seq, int layer, int index, WORD *buffer, const rectangle *area) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return read_opened_fits_partial(seq, layer, index, buffer, area);
		case SEQ_SER:
			return ser_read_opened_partial(seq->ser_file, layer, index, buffer, area);
		default:
			break;
	}
	return 0;
}


/*****************************************************************************
 *                         SEQUENCE DATA MANAGEMENT                          *
 * **************************************************************************/

/* if FWHM was calculated on the sequence, a minimisation exists for all
 * images, and when switching to a new image, it should be set as the only item
 * in the star list, in order to be displayed.
 * A special care is required in PSF_list.c:clear_stars_list(), to not free this data. */
void set_fwhm_star_as_star_list_with_layer(sequence *seq, int layer) {
	assert(seq->regparam);
	/* we chose here the first layer that has been allocated, which doesn't
	 * mean it contains data for all images. Handle with care. */
	if (seq->regparam && layer >= 0 && layer < seq->nb_layers && seq->regparam[layer] &&
			seq->regparam[layer][seq->current].fwhm_data && !com.stars) {
		com.stars = malloc(2 * sizeof(fitted_PSF *));
		com.stars[0] = seq->regparam[layer][seq->current].fwhm_data;
		com.stars[1] = NULL;
		com.star_is_seqdata = TRUE;
	}
}

// cannot be called in the worker thread
void set_fwhm_star_as_star_list(sequence *seq) {
	int layer = get_registration_layer();
	set_fwhm_star_as_star_list_with_layer(seq, layer);
}

/* Rebuilds the file name of an image in a sequence.
 * The file name is stored in name_buffer, which must be allocated 256 bytes
 * The index is the index in the sequence, not the number appearing in the file name
 * Return value: NULL on error, name_buffer on success.
 */
char *fit_sequence_get_image_filename(sequence *seq, int index, char *name_buffer, gboolean add_fits_ext) {
	char format[20];
	if (index < 0 || index > seq->number || name_buffer == NULL)
		return NULL;
	if (seq->fixed <= 1){
		sprintf(format, "%%s%%d");
	} else {
		sprintf(format, "%%s%%.%dd", seq->fixed);
	}
	if (add_fits_ext)
		strcat(format, com.ext);
	snprintf(name_buffer, 255, format,
			seq->seqname, seq->imgparam[index].filenum);
	name_buffer[255] = '\0';
	return name_buffer;
}

/* Returns a filename for an image that could be in a sequence, but the sequence structure
 * has not been fully initialized yet. Only beg, end, fixed and seqname are used.
 */
char *get_possible_image_filename(sequence *seq, int image_number, char *name_buffer) {
	char format[20];
	if (image_number < seq->beg || image_number > seq->end || name_buffer == NULL)
		return NULL;
	if (seq->fixed <= 1){
		sprintf(format, "%%s%%d%s", com.ext);
	} else {
		sprintf(format, "%%s%%.%dd%s", seq->fixed, com.ext);
	}
	sprintf(name_buffer, format, seq->seqname, image_number);
	return name_buffer;
}

/* splits a filename in a base name and an index number, if the file name ends with .fit
 * it also computes the fixed length if there are zeros in the index */
int	get_index_and_basename(const char *filename, char **basename, int *index, int *fixed){
	char *buffer;
	int i, fnlen, first_zero, digit_idx;

	*index = -1;		// error values
	*fixed = 0;
	first_zero = -1;
	*basename = NULL;
	fnlen = strlen(filename);
	if (fnlen < strlen(com.ext)+2) return -1;
	if (!ends_with(filename, com.ext)) return -1;
	i = fnlen-strlen(com.ext)-1;
	if (!isdigit(filename[i])) return -1;
	digit_idx = i;

	buffer = strdup(filename);
	buffer[fnlen-strlen(com.ext)] = '\0';		// for atoi()
	do {
		if (buffer[i] == '0' && first_zero < 0)
			first_zero = i;
		if (buffer[i] != '0' && first_zero > 0)
			first_zero = -1;
		i--;
	} while (i >= 0 && isdigit(buffer[i]));
	i++;
	if (i == 0) {
		free(buffer);
		return -1;	// no base name, only number
	}
	if (first_zero >= 0)
		*fixed = digit_idx - i + 1;
	//else *fixed = 0;
	*index = atoi(buffer+i);
	if (*basename == NULL) {	// don't copy it if we already have it
		*basename = malloc(i * sizeof(char) + 1);
		strncpy(*basename, buffer, i);
		(*basename)[i] = '\0';
	}
	//fprintf(stdout, "from filename %s, base name is %s, index is %d\n", filename, *basename, *index);
	free(buffer);
	return 0;
}

/* sets default values for the sequence */
void initialize_sequence(sequence *seq, gboolean is_zeroed) {
	int i;
	if (!is_zeroed) {
		memset(seq, 0, sizeof(sequence));
	}
	seq->nb_layers = -1;		// uninit value
	seq->reference_image = -1;	// uninit value
	seq->reference_star = -1;	// uninit value
	seq->type = SEQ_REGULAR;
	for (i=0; i<PREVIEW_NB; i++) {
		seq->previewX[i] = -1;
		seq->previewY[i] = -1;
	}
}

/* call this to close a sequence. Second arg must be FALSE for com.seq
 * WARNING: the data is not reset to NULL, if seq is to be reused,
 * initialize_sequence() must be called on it right after free_sequence()
 * (= do it for com.seq) */
void free_sequence(sequence *seq, gboolean free_seq_too) {
	static GtkComboBoxText *cbbt_layers = NULL;
	int i, j;
		
	if (cbbt_layers == NULL)
		cbbt_layers = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(
					builder, "comboboxreglayer"));
	gtk_combo_box_text_remove_all(cbbt_layers);
	
	if (seq == NULL) return;
	if (seq->nb_layers > 0 && seq->regparam) {
		for (i=0; i<seq->nb_layers; i++) {
			if (seq->regparam[i]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->regparam[i][j].fwhm_data
							&& ((seq->photometry[0] != NULL)
									&& seq->regparam[i][j].fwhm_data
											!= seq->photometry[0][j])) // avoid double free
						free(seq->regparam[i][j].fwhm_data);
				}
				free(seq->regparam[i]);
			}
		}
		free(seq->regparam);
	}

	for (i=0; i<seq->number; i++) {
		if (seq->fptr && seq->fptr[i]) {
			int status = 0;
			fits_close_file(seq->fptr[i], &status);
		}
		if (seq->imgparam) {
			if (seq->imgparam[i].stats)
				free(seq->imgparam[i].stats);
			if (seq->imgparam[i].date_obs)
				free(seq->imgparam[i].date_obs);
		}
	}
	if (seq->seqname)	free(seq->seqname);
	if (seq->layers)	free(seq->layers);
	if (seq->imgparam)	free(seq->imgparam);
	if (seq->fptr)		free(seq->fptr);

#ifdef _OPENMP
	if (seq->fd_lock) {
		for (i=0; i<seq->number; i++) {
			omp_destroy_lock(&seq->fd_lock[i]);
		}
		free(seq->fd_lock);
	}
#endif

		if (seq->ser_file) {
		ser_close_file(seq->ser_file);	// frees the data too
		free(seq->ser_file);
	}
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	if (seq->film_file) {
		film_close_file(seq->film_file);	// frees the data too
		free(seq->film_file);
	}
#endif

	if (seq->internal_fits) {
		/* the fits in internal_fits should still be referenced somewhere */
		free(seq->internal_fits);
	}
	/* Here this is a bit tricky. An internal sequence is a single image. So some
	 * processes like RGB alignment could free sequences and load it again: we need
	 * to keep undo history.
	 * In the case of a standard sequence, loading a new sequence MUST remove all
	 * undo history.
	 */
	if (seq->type != SEQ_INTERNAL)
		undo_flush();
	reset_plot();

	for (i = 0; i < MAX_SEQPSF && seq->photometry[i]; i++) {
		free_photometry_set(seq, i);
	}

	if (free_seq_too)	free(seq);
}

void sequence_free_preprocessing_data(sequence *seq) {
	// free opened files
	if (seq->ppprefix) {
		free(seq->ppprefix);
		seq->ppprefix = NULL;
	}
	if (seq->offset) {
		clearfits(seq->offset);
		free(seq->offset);
		seq->offset = NULL;
	}
	if (seq->dark) {
		clearfits(seq->dark);
		free(seq->dark);
		seq->dark = NULL;
	}
	if (seq->flat) {
		clearfits(seq->flat);
		free(seq->flat);
		seq->flat = NULL;
	}
}

gboolean sequence_is_loaded() {
	return (com.seq.seqname != NULL && com.seq.imgparam != NULL);
}

/* if no reference image has been set, return the index of an image that is
 * selected in the sequence, the best of the first registration data found if
 * any, the first otherwise */
int sequence_find_refimage(sequence *seq) {
	if (seq->reference_image != -1)
		return seq->reference_image;
	int layer, image, best = -1;
	for (layer = 0; layer < seq->nb_layers; layer++) {
		if (seq->regparam[layer]) {
			gboolean use_fwhm;
			double best_val;
			if (seq->regparam[layer][0].fwhm > 0.0) {
				use_fwhm = TRUE;
				best_val = 1000000.0;
			} else if (seq->regparam[layer][0].quality <= 0.0) {
				use_fwhm = FALSE;
				best_val = 0.0;
			}
			else continue;

			for (image = 0; image < seq->number; image++) {
				if (!seq->imgparam[image].incl)
					continue;
				if (use_fwhm) {
					if (seq->regparam[layer][image].fwhm > 0 &&
						       	seq->regparam[layer][image].fwhm < best_val) {
						best_val = seq->regparam[layer][image].fwhm;
						best = image;
					}
				} else {
					if (seq->regparam[layer][image].quality > 0 &&
						       	seq->regparam[layer][image].quality > best_val) {
						best_val = seq->regparam[layer][image].quality;
						best = image;
					}
				}
			}
		}
	}

	if (best == -1 && seq->selnum > 0) {
		for (image = 0; image < seq->number; image++) {
			if (seq->imgparam[image].incl) {
				best = image;
				break;
			}
		}
	}
	
	if (best == -1) best = 0;	// the first anyway if no regdata and no selected

	return best;
}

/* requires seq->nb_layers and seq->number to be already set */
void check_or_allocate_regparam(sequence *seq, int layer) {
	assert(layer < seq->nb_layers);
	if (!seq->regparam && seq->nb_layers > 0) {
		seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
		seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
	}
	if (seq->regparam && !seq->regparam[layer] && seq->number > 0) {
		seq->regparam[layer] = calloc(seq->number, sizeof(regdata));
	}
}

/* internal sequence are a set of 1-layer images already loaded elsewhere, and
 * directly referenced as fits *.
 * This is used in LRGV composition.
 * The returned sequence does not contain any reference to files, and thus has
 * to be populated with internal_sequence_set() */
sequence *create_internal_sequence(int size) {
	int i;
	sequence *seq = calloc(1, sizeof(sequence));
	initialize_sequence(seq, TRUE);
	seq->type = SEQ_INTERNAL;
	seq->number = size;
	seq->selnum = size;
	seq->nb_layers = 1;	
	seq->internal_fits = calloc(size, sizeof(fits *));
	seq->seqname = strdup(_("internal sequence"));
	seq->imgparam = calloc(size, sizeof(imgdata));
	for (i = 0; i < size; i++) {
		seq->imgparam[i].filenum = i;
		seq->imgparam[i].incl = 1;
		seq->imgparam[i].stats = NULL;
		seq->imgparam[i].date_obs = NULL;
}
check_or_allocate_regparam(seq, 0);
	return seq;
}

void internal_sequence_set(sequence *seq, int index, fits *fit) {
	assert(seq);
	assert(seq->internal_fits);
	assert(index < seq->number);
	seq->internal_fits[index] = fit;
}

// find index of the fit argument in the sequence
int internal_sequence_find_index(sequence *seq, fits *fit) {
	int i;
	assert(seq);
	assert(seq->internal_fits);
	for (i = 0; i < seq->number; i++) {
		if (fit == seq->internal_fits[i])
			return i;
	}
	return -1;
}

gboolean end_crop_sequence(gpointer p) {
	struct crop_sequence_data *args = (struct crop_sequence_data *) p;

	stop_processing_thread();// can it be done here in case there is no thread?
	if (!args->retvalue) {
		char *rseqname = malloc(
				strlen(args->prefix) + strlen(com.seq.seqname) + 5);

		sprintf(rseqname, "%s%s.seq", args->prefix, com.seq.seqname);
		check_seq(0);
		update_sequences_list(rseqname);
		free(rseqname);
	}
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

gpointer crop_sequence(gpointer p) {
	struct crop_sequence_data *args = (struct crop_sequence_data *) p;
	int frame, ret;
	float cur_nb;
	struct ser_struct *ser_file = NULL;

	args->retvalue = 0;

	if (args->seq->type == SEQ_SER) {
		char dest[256];

		ser_file = malloc(sizeof(struct ser_struct));
		if (ser_file == NULL) {
			printf("memory error: crop_sequence\n");
			args->retvalue = 1;
			gdk_threads_add_idle(end_crop_sequence, args);
		}
		sprintf(dest, "%s%s.ser", args->prefix, args->seq->seqname);
		if (ser_create_file(dest, ser_file, TRUE, com.seq.ser_file)) {
			siril_log_message(_("Creating the SER file failed, aborting.\n"));
			free(ser_file);
			ser_file = NULL;
			args->retvalue = 1;
			gdk_threads_add_idle(end_crop_sequence, args);
		}
	}

	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		if (!get_thread_run())
			break;
		ret = seq_read_frame(args->seq, frame, &(wfit[0]));
		if (!ret) {
			char dest[256], filename[256];

			crop(&(wfit[0]), args->area);
			switch (args->seq->type) {
			case SEQ_REGULAR:
				fit_sequence_get_image_filename(args->seq, frame, filename,
				TRUE);
				sprintf(dest, "%s%s", args->prefix, filename);
				savefits(dest, &wfit[0]);
				break;
			case SEQ_SER:
				if (ser_file) {
					ser_file->image_width = wfit[0].rx;
					ser_file->image_height = wfit[0].ry;
					if (ser_write_frame_from_fit(ser_file, &wfit[0], frame)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
					}
				}
				break;
			default:
				args->retvalue = 1;	// should not happen
			}

			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / args->seq->number);
		}
	}
	if (args->seq->type == SEQ_SER) {
		ser_write_and_close(ser_file);
		free(ser_file);
	}
	gdk_threads_add_idle(end_crop_sequence, args);
	return 0;
}

// check if the passed sequence is used as a color sequence. It can be a CFA
// sequence explicitly demoisaiced too, which returns true.
gboolean sequence_is_rgb(sequence *seq) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return seq->nb_layers == 3;
		case SEQ_SER:
			return (seq->ser_file->color_id != SER_MONO && com.debayer.open_debayer) ||
				seq->ser_file->color_id == SER_RGB ||
				seq->ser_file->color_id == SER_BGR;
		default:
			return TRUE;
	}
}

/* Get statistics for an image in a sequence.
 * If it's not in the cache, it will be computed from the_image. If the_image is NULL,
 * it returns NULL in that case.
 * Do not free result.
 */
imstats* seq_get_imstats(sequence *seq, int index, fits *the_image, int option) {
	assert(seq->imgparam);
	if (!seq->imgparam[index].stats && the_image) {
		seq->imgparam[index].stats = statistics(the_image, 0, NULL, option, STATS_ZERO_NULLCHECK);
		if (!seq->imgparam[index].stats) {
			siril_log_message(_("Error: no data computed.\n"));
			return NULL;
		}
		seq->needs_saving = TRUE;
	}
	return seq->imgparam[index].stats;
}

/* Ensures that an area does not derive off-image.
 * Verifies coordinates of the center and moves it inside the image if the area crosses the bounds.
 */
void enforce_area_in_image(rectangle *area, sequence *seq) {
	if (area->x < 0) area->x = 0;
	if (area->y < 0) area->y = 0;
	if (area->x + area->w > seq->rx)
		area->x = seq->rx - area->w;
	if (area->y + area->h > seq->ry)
		area->y = seq->ry - area->h;
}

struct exportseq_args {
	sequence *seq;
	char *basename;
	int convflags;
	gboolean normalize;
//	int gif_delay, gif_loops;
	double avi_fps;
	int quality;	// [1, 5], for mp4 and webm
	gboolean resize;
	int32_t dest_width, dest_height;
	gboolean crop;
	rectangle crop_area;
	seq_image_filter filtering_criterion;
	double filtering_parameter;
};

/* Used for avi exporter */
static uint8_t *fits_to_uint8(fits *fit) {
	uint8_t *data;
	int w, h, i, j, channel, step;
	float pente;
	WORD lo, hi;

	w = fit->rx;
	h = fit->ry;
	channel = fit->naxes[2];
	step = (channel == 3 ? 2 : 0);
	pente = computePente(&lo, &hi);

	data = malloc(w * h * channel * sizeof(uint8_t));
	for (i = 0, j = 0; i < w * h * channel; i += channel, j++) {
		data[i + step] = (uint8_t) round_to_BYTE(((double) fit->pdata[RLAYER][j] * pente));
		if (channel > 1) {
			data[i + 1] = (uint8_t) round_to_BYTE(((double) fit->pdata[GLAYER][j] * pente));
			data[i + 2 - step] = (uint8_t) round_to_BYTE(((double) fit->pdata[BLAYER][j] * pente));
		}
	}
	return data;
}

gpointer export_sequence(gpointer ptr) {
	int i, x, y, nx, ny, shiftx, shifty, layer, retval = 0, reglayer, nb_layers, skipped;
	float cur_nb = 0.f, nb_frames;
	unsigned int out_width, out_height, in_width, in_height, nbdata = 0;
	uint8_t *data;
	fits fit, destfit;
	char filename[256], dest[256];
	struct ser_struct *ser_file = NULL;
	GSList *timestamp = NULL;
	gchar *strTime;
#ifdef HAVE_FFMPEG
	struct mp4_struct *mp4_file = NULL;
#endif
	struct exportseq_args *args = (struct exportseq_args *)ptr;
	norm_coeff coeff;
	memset(&fit, 0, sizeof(fits));
	memset(&destfit, 0, sizeof(fits));

	reglayer = get_registration_layer();
	siril_log_message(_("Using registration information from layer %d to export sequence\n"), reglayer);
	if (args->crop) {
		in_width  = args->crop_area.w;
		in_height = args->crop_area.h;
	} else {
		in_width  = args->seq->rx;
		in_height = args->seq->ry;
	}

	if (args->resize) {
		out_width = args->dest_width;
		out_height = args->dest_height;
		if (out_width == in_width && out_height == in_height)
			args->resize = FALSE;
	} else {
		out_width = in_width;
		out_height = in_height;
	}

	switch (args->convflags) {
		case TYPESER:
			/* image size is not known here, no problem for crop or resize */
			ser_file = malloc(sizeof(struct ser_struct));
			snprintf(dest, 256, "%s.ser", args->basename);
			if (ser_create_file(dest, ser_file, TRUE, NULL))
				siril_log_message(_("Creating the SER file failed, aborting.\n"));
			break;

		case TYPEAVI:
			/* image size is set here, resize is managed by opencv when
			 * writing frames, we don't need crop size here */
			snprintf(dest, 256, "%s.avi", args->basename);
			int32_t avi_format;

			if (args->seq->nb_layers == 1)
				avi_format = AVI_WRITER_INPUT_FORMAT_MONOCHROME;
			else avi_format = AVI_WRITER_INPUT_FORMAT_COLOUR;

			avi_file_create(dest, out_width, out_height, avi_format,
					AVI_WRITER_CODEC_DIB, args->avi_fps);
			break;

		case TYPEMP4:
		case TYPEWEBM:
#ifndef HAVE_FFMPEG
			siril_log_message(_("MP4 output is not supported because siril was not compiled with ffmpeg support.\n"));
			retval = -1;
			goto free_and_reset_progress_bar;
#else
			/* image size is set here, resize is managed by ffmpeg so it also
			 * needs to know the input image size after crop */
			snprintf(dest, 256, "%s.%s", args->basename,
					args->convflags == TYPEMP4 ? "mp4" : "webm");
			if (args->avi_fps <= 0) args->avi_fps = 25;

			if (in_width % 32 || out_height % 2 || out_width % 2) {
				siril_log_message(_("Film output needs to have a width that is a multiple of 32 and an even height, resizing selection.\n"));
				if (in_width % 32) in_width = (in_width / 32) * 32 + 32;
				if (in_height % 2) in_height++;
				if (args->crop) {
					args->crop_area.w = in_width;
					args->crop_area.h = in_height;
				} else {
					args->crop = TRUE;
					args->crop_area.x = 0;
					args->crop_area.y = 0;
					args->crop_area.w = in_width;
					args->crop_area.h = in_height;
				}
				compute_fitting_selection(&args->crop_area, 32, 2, 0);
				memcpy(&com.selection, &args->crop_area, sizeof(rectangle));
				fprintf(stdout, "final input area: %d,%d,\t%dx%d\n",
						args->crop_area.x, args->crop_area.y,
						args->crop_area.w, args->crop_area.h);
				in_width = args->crop_area.w;
				in_height = args->crop_area.h;
				if (!args->resize) {
					out_width = in_width;
					out_height = in_height;
				} else {
					if (out_width % 2) out_width++;
					if (out_height % 2) out_height++;
				}
			}

			mp4_file = mp4_create(dest, out_width, out_height, args->avi_fps, args->seq->nb_layers, args->quality, in_width, in_height);
			if (!mp4_file) {
				retval = -1;
				goto free_and_reset_progress_bar;
			}
#endif
			break;
	}

	if (args->normalize) {
		struct stacking_args stackargs;

		coeff.offset = malloc(args->seq->number * sizeof(double));
		// mul is not used in ADDITIVE_SCALING but needed to avoid crash in compute_normalization
		coeff.mul = malloc(args->seq->number * sizeof(double));
		coeff.scale = malloc(args->seq->number * sizeof(double));

		stackargs.force_norm = FALSE;
		stackargs.seq = args->seq;
		stackargs.nb_images_to_stack = args->seq->number;
		stackargs.filtering_criterion = stack_filter_all;

		stackargs.image_indices = malloc(stackargs.nb_images_to_stack * sizeof(int));
		fill_list_of_unfiltered_images(&stackargs);
		compute_normalization(&stackargs, &coeff, ADDITIVE_SCALING);
		// the image_indices are not used in the rest of this function for now
		free(stackargs.image_indices);
		if (args->seq->needs_saving)	// if we had to compute new stats
			writeseqfile(args->seq);
	}

	nb_frames = (float)args->seq->number;

	set_progress_bar_data(NULL, PROGRESS_RESET);
	for (i = 0, skipped = 0; i < args->seq->number; ++i) {
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, i, args->filtering_parameter)) {
			siril_log_message(_("image %d is excluded from export\n"), i);
			skipped++;
			continue;
		}

		if (!seq_get_image_filename(args->seq, i, filename)) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		char *tmpmsg = strdup(_("Processing image "));
		tmpmsg = str_append(&tmpmsg, filename);
		set_progress_bar_data(tmpmsg,
				(double) cur_nb / ((double) nb_frames + 1.));
		free(tmpmsg);

		if (seq_read_frame(args->seq, i, &fit)) {
			siril_log_message(_("Export: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		if (!nbdata) {
			/* destfit is allocated to the real size because of the possible
			 * shifts and of the inplace cropping. FITS data is copied from
			 * fit, image buffers are duplicated. */
			memcpy(&destfit, &fit, sizeof(fits));
			destfit.header = NULL;
			destfit.fptr = NULL;
			nbdata = fit.rx * fit.ry;
			destfit.data = calloc(nbdata * fit.naxes[2], sizeof(WORD));
			if (!destfit.data) {
				fprintf(stderr, "Could not allocate memory for the export, aborting\n");
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			destfit.pdata[0] = destfit.data;
			if (fit.naxes[2] == 1) {
				destfit.pdata[1] = destfit.data;
				destfit.pdata[2] = destfit.data;
			} else {
				destfit.pdata[1] = destfit.data + nbdata;
				destfit.pdata[2] = destfit.data + nbdata * 2;
			}
			nb_layers = fit.naxes[2];
		}
		else if (fit.ry * fit.rx != nbdata || nb_layers != fit.naxes[2]) {
			fprintf(stderr, "An image of the sequence doesn't have the same dimensions\n");
			retval = -3;
			goto free_and_reset_progress_bar;
		}
		else {
			/* we want copy the header */
			copy_header(&fit, &destfit);
			memset(destfit.data, 0, nbdata * fit.naxes[2] * sizeof(WORD));
			if (args->crop) {
				/* reset destfit damaged by the crop function */
				if (fit.naxes[2] == 3) {
					destfit.pdata[1] = destfit.data + nbdata;
					destfit.pdata[2] = destfit.data + nbdata * 2;
				}
				destfit.rx = destfit.naxes[0] = fit.rx;
				destfit.ry = destfit.naxes[1] = fit.ry;
			}
		}

		/* load registration data for current image */
		if (reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = args->seq->regparam[reglayer][i].shiftx;
			shifty = args->seq->regparam[reglayer][i].shifty;
		} else {
			shiftx = 0;
			shifty = 0;
		}

		/* fill the image with shift data and normalization */
		for (layer=0; layer<fit.naxes[2]; ++layer) {
			for (y=0; y < fit.ry; ++y){
				for (x=0; x < fit.rx; ++x){
					nx = x + shiftx;
					ny = y + shifty;
					if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
						if (args->normalize) {
							double tmp = fit.pdata[layer][x + y * fit.rx];
							tmp *= coeff.scale[i];
							tmp -= coeff.offset[i];
							destfit.pdata[layer][nx + ny * fit.rx] = round_to_WORD(tmp);
						} else {
							destfit.pdata[layer][nx + ny * fit.rx] = fit.pdata[layer][x + y * fit.rx];
						}
					}
				}
			}
		}

		if (args->crop) {
			crop(&destfit, &args->crop_area);
		}

		switch (args->convflags) {
			case TYPEFITS:
				snprintf(dest, 255, "%s%05d%s", args->basename, i, com.ext);
				if (savefits(dest, &destfit)) {
					retval = -1;
					goto free_and_reset_progress_bar;
				}
				break;
			case TYPESER:
				strTime = strdup(destfit.date_obs);
				timestamp = g_slist_append (timestamp, strTime);
				if (ser_write_frame_from_fit(ser_file, &destfit, i - skipped))
					siril_log_message(
							_("Error while converting to SER (no space left?)\n"));
				break;
			case TYPEAVI:
				data = fits_to_uint8(&destfit);

				if (args->resize) {
#ifdef HAVE_OPENCV
					uint8_t *newdata = malloc(out_width * out_height * destfit.naxes[2]);
					cvResizeGaussian_data8(data, destfit.rx, destfit.ry, newdata,
							out_width, out_height, destfit.naxes[2], OPENCV_CUBIC);
					avi_file_write_frame(0, newdata);
					free(newdata);
#else
					siril_log_message(_("Siril needs opencv to resize images\n"));
					avi_file_write_frame(0, data);
#endif
				}
				else
					avi_file_write_frame(0, data);
				free(data);
				break;
#ifdef HAVE_FFMPEG
			case TYPEMP4:
			case TYPEWEBM:
				mp4_add_frame(mp4_file, &destfit);
				break;
#endif
		}
		cur_nb += 1.f;
		set_progress_bar_data(NULL, cur_nb / nb_frames);

		clearfits(&fit);
	}

free_and_reset_progress_bar:
	clearfits(&fit);	// in case of goto
	clearfits(&destfit);
	if (args->normalize) {
		free(coeff.offset);
		free(coeff.mul);
		free(coeff.scale);
	}
	if (args->convflags == TYPESER) {
		ser_convertTimeStamp(ser_file, timestamp);
		ser_write_and_close(ser_file);
		free(ser_file);
		g_slist_free_full(timestamp, g_free);
	}
	else if (args->convflags == TYPEAVI) {
		avi_file_close(0);
	}
#ifdef HAVE_FFMPEG
	else if (mp4_file && (args->convflags == TYPEMP4 || args->convflags == TYPEWEBM)) {
		mp4_close(mp4_file);
		free(mp4_file);
	}
#endif

	if (retval) {
		set_progress_bar_data(_("Sequence export failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Sequence export failed\n"));
	}
	else {
		set_progress_bar_data(_("Sequence export succeeded."), PROGRESS_RESET);
		siril_log_message(_("Sequence export succeeded.\n"));
	}

	free(args->basename);
	free(args);
	args = NULL;
	gdk_threads_add_idle(end_generic, args);
	return NULL;
}

void on_buttonExportSeq_clicked(GtkButton *button, gpointer user_data) {
	int selected = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboExport")));
	const char *bname = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryExportSeq")));
	GtkAdjustment *stackadj = GTK_ADJUSTMENT(gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(lookup_widget("stackspin"))));
	GtkComboBox *stack_type = GTK_COMBO_BOX(lookup_widget("comboboxstacksel"));
	struct exportseq_args *args;
	GtkToggleButton *exportNormalize, *checkResize;
	GtkEntry *fpsEntry, *widthEntry, *heightEntry;
	GtkAdjustment *adjQual;
	double percent;

	if (bname[0] == '\0') return;
	if (selected == -1) return;

	args = malloc(sizeof(struct exportseq_args));
	args->basename = strdup(bname);
	args->seq = &com.seq;
	exportNormalize = GTK_TOGGLE_BUTTON(lookup_widget("exportNormalize"));
	args->normalize = gtk_toggle_button_get_active(exportNormalize);
	args->crop = com.selection.w && com.selection.h;
	if (args->crop)
		memcpy(&args->crop_area, &com.selection, sizeof(rectangle));

	// filtering
	switch (gtk_combo_box_get_active(stack_type)) {
		case 0:
			args->filtering_criterion = stack_filter_all;
			break;
		case 1:
			args->filtering_criterion = stack_filter_included;
			break;
		case 2:
			percent = gtk_adjustment_get_value(stackadj);
			args->filtering_criterion = stack_filter_fwhm;
			args->filtering_parameter = compute_highest_accepted_fwhm(percent);
			break;
		case 3:
			percent = gtk_adjustment_get_value(stackadj);
			args->filtering_criterion = stack_filter_quality;
			args->filtering_parameter = compute_highest_accepted_quality(percent);
	}

	// format
	switch (selected) {
	case 0:
		args->convflags = TYPEFITS;
		args->basename = format_basename(args->basename);
		break;
	case 1:
		args->convflags = TYPESER;
		break;
	case 2:
	case 3:
	case 4:
		fpsEntry = GTK_ENTRY(lookup_widget("entryAviFps"));
		args->avi_fps = atoi(gtk_entry_get_text(fpsEntry));
		widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));
		args->dest_width = atof(gtk_entry_get_text(widthEntry));
		heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));
		args->dest_height = atof(gtk_entry_get_text(heightEntry));
		checkResize = GTK_TOGGLE_BUTTON(lookup_widget("checkAviResize"));
		adjQual = GTK_ADJUSTMENT(gtk_builder_get_object(builder,"adjustment3"));
		args->quality = (int)gtk_adjustment_get_value(adjQual);

		if (args->dest_height == 0 || args->dest_width == 0) {
			siril_log_message(_("Width or height cannot be null. Not resizing.\n"));
			args->resize = FALSE;
			gtk_toggle_button_set_active(checkResize, FALSE);
		} else if (args->dest_height == args->seq->ry && args->dest_width == args->seq->rx) {
			args->resize = FALSE;
			gtk_toggle_button_set_active(checkResize, FALSE);
		} else {
			args->resize = gtk_toggle_button_get_active(checkResize);
		}
		args->convflags = TYPEAVI;
		if (selected == 3)
			args->convflags = TYPEMP4;
		else if (selected == 4)
			args->convflags = TYPEWEBM;
		break;
	default:
		free(args);
		return;
	}
	set_cursor_waiting(TRUE);
	start_in_new_thread(export_sequence, args);
}

void on_comboExport_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *avi_options = lookup_widget("boxAviOptions");
	GtkWidget *checkAviResize = lookup_widget("checkAviResize");
	GtkWidget *quality = lookup_widget("exportQualScale");
	gtk_widget_set_visible(avi_options, gtk_combo_box_get_active(box) >= 2);
	gtk_widget_set_visible(quality, gtk_combo_box_get_active(box) >= 3);
#ifdef HAVE_OPENCV
	gtk_widget_set_sensitive(checkAviResize, TRUE);
#else
	gtk_widget_set_sensitive(checkAviResize, FALSE);
#endif
}

void on_checkAviResize_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *heightEntry = lookup_widget("entryAviHeight");
	GtkWidget *widthEntry = lookup_widget("entryAviWidth");
	gtk_widget_set_sensitive(heightEntry, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(widthEntry, gtk_toggle_button_get_active(togglebutton));
}

void update_export_crop_label() {
	static GtkLabel *label = NULL;
	if (!label) 
		label = GTK_LABEL(lookup_widget("exportLabel"));
	if (com.selection.w && com.selection.h)
		gtk_label_set_text(label, _("Cropping to selection"));
	else gtk_label_set_text(label, _("Select area to crop"));
}

/********************************************************************
 *                                             __                   *
 *                   ___  ___  __ _ _ __  ___ / _|                  *
 *                  / __|/ _ \/ _` | '_ \/ __| |_                   *
 *                  \__ \  __/ (_| | |_) \__ \  _|                  *
 *                  |___/\___|\__, | .__/|___/_|                    *
 *                               |_|_|                              *
 ********************************************************************/

struct seqpsf_args {
	gboolean for_registration;
	framing_mode framing;

	/* The seqpsf result for each image */
	GSList *list;
};

struct seqpsf_data {
	int image_index;
	fitted_PSF *psf;
	double exposure;
};

/* Computes FWHM for a sequence image.
 * area is the area from which fit was extracted from the full frame.
 * when the framing is set to follow star, args->area is centered on the found star
 */
int seqpsf_image_hook(struct generic_seq_args *args, int index, fits *fit, rectangle *area) {
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;
	struct seqpsf_data *data = malloc(sizeof(struct seqpsf_data));
	data->image_index = index;

	rectangle psfarea = { .x = 0, .y = 0, .w = fit->rx, .h = fit->ry };
	data->psf = psf_get_minimisation(fit, 0, &psfarea, !spsfargs->for_registration);
	if (data->psf) {
		data->psf->xpos = data->psf->x0 + area->x;
		data->psf->ypos = area->y + area->h - data->psf->y0;

		/* let's move args->area to center it on the star */
		if (spsfargs->framing == FOLLOW_STAR_FRAME) {
			args->area.x = round_to_int(data->psf->xpos) - args->area.w/2;
			args->area.y = round_to_int(data->psf->ypos) - args->area.h/2;
		}

		if (!args->seq->imgparam[index].date_obs && fit->date_obs[0] != '\0')
			args->seq->imgparam[index].date_obs = strdup(fit->date_obs);
		data->exposure = fit->exposure;
	}
	else {
		siril_log_color_message(_("No star found in the area image %d around %d,%d\n"),
				"red", index, area->x, area->y);
	}

#ifdef _OPENMP
	omp_set_lock(&args->lock);
#endif
	spsfargs->list = g_slist_prepend(spsfargs->list, data);
#ifdef _OPENMP
	omp_unset_lock(&args->lock);
#endif
	return !data->psf;
}

gboolean end_seqpsf(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *)p;
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;
	sequence *seq = args->seq;
	int layer = args->layer_for_partial;
	int photometry_index;
	gboolean displayed_warning = FALSE, write_to_regdata = FALSE;
	gboolean dont_stop_thread;

	if (args->retval)
		goto proper_ending;

	if (spsfargs->for_registration || !seq->regparam[layer]) {
		check_or_allocate_regparam(seq, layer);
		write_to_regdata = TRUE;
	}
	if (!spsfargs->for_registration) {
		int i;
		for (i = 0; i < MAX_SEQPSF && seq->photometry[i]; i++);
		if (i == MAX_SEQPSF) {
			free_photometry_set(seq, 0);
		       	i = 0;
		}
		seq->photometry[i] = calloc(seq->number, sizeof(fitted_PSF *));
		photometry_index = i;
	}

	GSList *iterator = spsfargs->list;
	while (iterator) {
		struct seqpsf_data *data = iterator->data;
		/* check exposure consistency (only obtained when !for_registration) */
		if (!spsfargs->for_registration && seq->exposure > 0.0 &&
				seq->exposure != data->exposure && !displayed_warning) {
			siril_log_color_message(_("Star analysis does not give consistent results when exposure changes across the sequence.\n"), "red");
			displayed_warning = TRUE;
		}
		seq->exposure = data->exposure;

		if (write_to_regdata) {
			seq->regparam[layer][data->image_index].fwhm_data = data->psf;
			if (data->psf)
				seq->regparam[layer][data->image_index].fwhm = data->psf->fwhmx;
		}

		// for photometry use: store data in seq->photometry
		if (!spsfargs->for_registration) {
			seq->photometry[photometry_index][data->image_index] = data->psf;
		}

		iterator = g_slist_next(iterator);
	}

	// for registration use: store data in seq->regparam
	if (spsfargs->for_registration || !seq->regparam[layer]) {
		// the data put in regparam if !for_registration is not yet saved
		seq->needs_saving = TRUE;
	}

	if (com.seq.needs_saving)
		writeseqfile(&com.seq);

	set_fwhm_star_as_star_list_with_layer(seq, layer);

	if (!args->already_in_a_thread) {
		/* do here all GUI-related items, because it runs in the main thread.
		 * Most of these things are already done in end_register_idle
		 * in case seqpsf is called for registration. */
		// update the list in the GUI
		if (seq->type != SEQ_INTERNAL)
			fill_sequence_list(seq, layer);

		set_layers_for_registration();	// update display of available reg data
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
	}

proper_ending:
	dont_stop_thread = args->already_in_a_thread;
	if (spsfargs->list)
		g_slist_free(spsfargs->list);
	free(spsfargs);

	if (dont_stop_thread) {
		// we must not call stop_processing_thread() here
		return FALSE;
	} else {
		free(args);
		return end_generic(NULL);
	}
}

/* process PSF for the given sequence, on the given layer, the area of the
 * image selection (com.selection), as a threaded operation or not.
 */
int seqpsf(sequence *seq, int layer, gboolean for_registration,
		framing_mode framing, gboolean run_in_thread) {

	if (framing == REGISTERED_FRAME && !seq->regparam[layer])
		framing = ORIGINAL_FRAME;

	if (com.selection.w <= 0 || com.selection.h <= 0){
		siril_log_message(_("Select an area first\n"));
		return 1;
	}

	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	struct seqpsf_args *spsfargs = malloc(sizeof(struct seqpsf_args));

	spsfargs->for_registration = for_registration;
	spsfargs->framing = framing;
	spsfargs->list = NULL;	// GSList init is NULL

	args->seq = seq;
	args->partial_image = TRUE;
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	args->layer_for_partial = layer;
	args->regdata_for_partial = framing == REGISTERED_FRAME;
	args->get_photometry_data_for_partial = !for_registration;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = seq->selnum;
	args->prepare_hook = NULL;
	args->finalize_hook = NULL;
	args->image_hook = seqpsf_image_hook;
	args->idle_function = end_seqpsf;
	args->description = _("PSF on area");
	args->has_output = FALSE;
	args->user = spsfargs;
	args->already_in_a_thread = !run_in_thread;
	args->parallel = framing != FOLLOW_STAR_FRAME;

	if (run_in_thread) {
		start_in_new_thread(generic_sequence_worker, args);
		return 0;
	} else {
		generic_sequence_worker(args);
		int retval = args->retval;
		free(args);
		return retval;
	}
}

