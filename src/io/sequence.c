/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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

#ifdef _WIN32
#include <windows.h>
#endif

#include <gtk/gtk.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <libgen.h>

#include "core/siril.h"
#include "sequence.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/plot.h"
#include "ser.h"
#ifdef HAVE_FFMS2
#include "films.h"
#endif
#include "avi_pipp/avi_writer.h"
#include "single_image.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"	// clear_stars_list
#include "algos/PSF.h"
#include "algos/quality.h"
#include "algos/statistics.h"
#include "algos/geometry.h"
#include "registration/registration.h"
#include "stacking/stacking.h"	// for update_stack_interface


/* com.seq is a static struct containing the sequence currently selected by the
 * user from the interface. It may change to be a pointer to any sequence
 * someday, until then, the seqname is NULL when no sequence is loaded and the
 * number of images in the sequence is also negative.
 * com.uniq represents information about an image opened and displayed outside
 * a sequence, for example from the load command, the open menu, or the result
 * of a stacking operation.
 * com.seq.number is used to provide a relationship between a possibly loaded
 * sequence and the single image. A single image can be loaded without
 * unloading the sequence. This information could as well be moved to
 * com.status if com.seq becomes a pointer. Three constants have been declared
 * in siril.h to explicit this relationship: RESULT_IMAGE, UNRELATED_IMAGE and
 * SCALED_IMAGE. They are mostly used to understand what to do to display
 * single images when a sequence is loaded or not.
 */

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
	if (!changedir(dirname, NULL)) {
		writeinitfile();
		if (!com.script) {
			set_GUI_CWD();
		}
	}
	g_free(dirname);

	if (check_only_one_film_seq(realname)) retval = 1;
	else {
#ifdef HAVE_FFMS2
		const char *ext;
#endif
		switch (imagetype) {
			case TYPESER:
				name[strlen(name)-1] = 'q';
				break;
#ifdef HAVE_FFMS2
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
 * force clears the stats in the seqfile.
 */
int check_seq(int recompute_stats) {
	char *basename;
	int curidx, fixed;
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	sequence **sequences;
	int i, nb_seq = 0, max_seq = 10;

	if (!com.wd) {
		siril_log_message(_("Current working directory is not set, aborting.\n"));
		return 1;
	}
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		fprintf (stderr, "check_seq: %s\n", error->message);
		g_error_free(error);
		g_free(com.wd);
		com.wd = NULL;
		return 1;
	}

	sequences = malloc(sizeof(sequence *) * max_seq);
	set_progress_bar_data(NULL, PROGRESS_PULSATE);

	while ((file = g_dir_read_name(dir)) != NULL) {
		sequence *new_seq;
		int fnlen = strlen(file);
		if (fnlen < 4) continue;
		const char *ext = get_filename_ext(file);
		if (!ext) continue;
		if (!strcasecmp(ext, "ser")) {
			struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
			ser_init_struct(ser_file);
			if (ser_open_file(file, ser_file))
				continue;
			new_seq = calloc(1, sizeof(sequence));
			initialize_sequence(new_seq, TRUE);
			new_seq->seqname = g_strndup(file, fnlen-4);
			new_seq->beg = 0;
			new_seq->end = ser_file->frame_count-1;
			new_seq->number = ser_file->frame_count;
			new_seq->type = SEQ_SER;
			new_seq->ser_file = ser_file;
			sequences[nb_seq] = new_seq;
			nb_seq++;
			fprintf(stdout, "Found a SER sequence (number %d)\n", nb_seq);
		}
#ifdef HAVE_FFMS2
		else if (!check_for_film_extensions(ext)) {
			struct film_struct *film_file = malloc(sizeof(struct film_struct));
			if (film_open_file(file, film_file)) {
				free(film_file);
				continue;
			}
			new_seq = calloc(1, sizeof(sequence));
			initialize_sequence(new_seq, TRUE);
			int len = strlen(ext);
			new_seq->seqname = g_strndup(file, fnlen - (len + 1));
			new_seq->beg = 0;
			new_seq->end = film_file->frame_count-1;
			new_seq->number = film_file->frame_count;
			new_seq->type = SEQ_AVI;
			new_seq->film_file = film_file;
			sequences[nb_seq] = new_seq;
			nb_seq++;
			fprintf(stdout, "Found a AVI sequence (number %d)\n", nb_seq);
		}
#endif

		else if (!strcasecmp(ext, com.ext+1)) {
			if (!get_index_and_basename(file, &basename, &curidx, &fixed)) {
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
				PRINT_ALLOC_ERR;
				break;
			}
		}
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);
	g_dir_close(dir);
	if (nb_seq > 0) {
		int retval = 1;
		for (i=0; i<nb_seq; i++) {
			if (sequences[i]->beg != sequences[i]->end) {
				char msg[200];
				sprintf(msg, _("sequence %d, found: %d to %d"),
						i+1, sequences[i]->beg, sequences[i]->end);
				if (!buildseqfile(sequences[i], recompute_stats) && retval)
					retval = 0;	// at least one succeeded to be created
			}
			free_sequence(sequences[i], TRUE);
		}
		free(sequences);
		return retval;
	}
	free(sequences);
	siril_log_message(_("No sequence found, verify working directory or "
			"change FITS extension in settings (current is %s)\n"), com.ext);
	return 1;	// no sequence found
}

/* Check for on film sequence of the name passed in arguement
 * Returns 0 if OK */
int check_only_one_film_seq(char* name) {
	int retval = 1;
	GDir *dir;
	GError *error = NULL;
	sequence *new_seq = NULL;

	if (!com.wd) {
		siril_log_message(_("Current working directory is not set, aborting.\n"));
		return 1;
	}
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		fprintf (stderr, "check_only_one_film_seq: %s\n", error->message);
		g_error_free(error);
		g_free(com.wd);
		com.wd = NULL;
		return 1;
	}

	int fnlen = strlen(name);
	const char *ext = get_filename_ext(name);

	if (!strcasecmp(ext, "ser")) {
		struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
		ser_init_struct(ser_file);
		if (ser_open_file(name, ser_file)) {
			g_dir_close(dir);
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
#ifdef HAVE_FFMS2
	else if (!check_for_film_extensions(ext)) {
		struct film_struct *film_file = malloc(sizeof(struct film_struct));
		if (film_open_file(name, film_file)) {
			free(film_file);
			g_dir_close(dir);
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
	g_dir_close(dir);
	if (!new_seq) return 0;
	if (new_seq->beg != new_seq->end) {
		if (!buildseqfile(new_seq, 0) && retval)
			retval = 0;
	}
	free_sequence(new_seq, TRUE);
	return retval;
}

// get the number of layers and image size for a new sequence
// if load_ref_into_gfit is true, the image is kept into gfit if data loading was
// required, and 1 is returned when it required loading
int seq_check_basic_data(sequence *seq, gboolean load_ref_into_gfit) {
	if (seq->nb_layers == -1 || seq->rx == 0) {	// not init yet, first loading of the sequence
		int image_to_load = sequence_find_refimage(seq);
		fits tmpfit = { 0 }, *fit;

		if (load_ref_into_gfit) {
			clearfits(&gfit);
			fit = &gfit;
		} else {
			fit = &tmpfit;
			memset(fit, 0, sizeof(fits));
		}

		if (seq_read_frame(seq, image_to_load, fit)) {
			fprintf(stderr, "could not load first image from sequence\n");
			return -1;
		}

		/* initialize sequence-related runtime data */
		seq->rx = fit->rx; seq->ry = fit->ry;
		seq->bitpix = fit->orig_bitpix;	// for partial read
		seq->data_max = fit->data_max; // for partial read
		fprintf(stdout, "bitpix for the sequence is set as %d\n", seq->bitpix);
		if (seq->nb_layers == -1) {
			seq->nb_layers = fit->naxes[2];
			seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
			seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
		}
		seq->needs_saving = TRUE;

		if (load_ref_into_gfit) {
			seq->current = image_to_load;
		} else {
			clearfits(fit);
		}
		return 1;
	}
	return 0;
}

static void free_cbbt_layers() {
	GtkComboBoxText *cbbt_layers = GTK_COMBO_BOX_TEXT(lookup_widget("comboboxreglayer"));
	gtk_combo_box_text_remove_all(cbbt_layers);
}

/* load a sequence and initializes everything that relates */
int set_seq(const char *name){
	sequence *seq = NULL;
	char *basename;

	if ((seq = readseqfile(name)) == NULL) {
		fprintf(stderr, "could not load sequence %s\n", name);
		return 1;
	}
	free_image_data();

	int retval = seq_check_basic_data(seq, TRUE);
	if (retval == -1) {
		free(seq);
		return 1;
	}
	if (retval == 0) {
		int image_to_load = sequence_find_refimage(seq);
		if (seq_read_frame(seq, image_to_load, &gfit)) {
			fprintf(stderr, "could not load first image from sequence\n");
			free(seq);
			return 1;
		}
		seq->current = image_to_load;
	}

	basename = g_path_get_basename(seq->seqname);
	siril_log_message(_("Sequence loaded: %s (%d->%d)\n"),
			basename, seq->beg, seq->end);
	g_free(basename);

	/* Sequence is stored in com.seq for now */
	close_sequence(TRUE);
	memcpy(&com.seq, seq, sizeof(sequence));

	if (seq->nb_layers > 1)
		show_rgb_window();
	else hide_rgb_window();
	init_layers_hi_and_lo_values(MIPSLOHI); // set some hi and lo values in seq->layers,
	set_cutoff_sliders_max_values();// update min and max values for contrast sliders
	set_cutoff_sliders_values();	// update values for contrast sliders for this image
	seqsetnum(seq->current);	// set limits for spin button and display loaded filenum
	set_layers_for_assign();	// set default layers assign and populate combo box
	set_layers_for_registration();	// set layers in the combo box for registration
	fill_sequence_list(seq, 0, FALSE);// display list of files in the sequence
	set_output_filename_to_sequence_name();
	sliders_mode_set_state(com.sliders);
	initialize_display_mode();
	reset_plot(); // reset all plots

	/* initialize image-related runtime data */
	set_display_mode();		// display the display mode in the combo box
	display_filename();		// display filename in gray window
	adjust_exclude(seq->current, FALSE);	// check or uncheck excluded checkbox
	adjust_refimage(seq->current);	// check or uncheck reference image checkbox
	update_prepro_interface(seq->type == SEQ_REGULAR); // enable or not the preprobutton
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
int seq_load_image(sequence *seq, int index, gboolean load_it) {
	if (!single_image_is_loaded())
		save_stats_from_fit(&gfit, seq, seq->current);
	clear_stars_list();
	clear_histograms();
	undo_flush();
	close_single_image();
	clearfits(&gfit);
	if (seq->current == SCALED_IMAGE) {
		gfit.rx = seq->rx;
		gfit.ry = seq->ry;
		adjust_vport_size_to_image();
	}
	seq->current = index;

	if (load_it) {
		set_cursor_waiting(TRUE);
		if (seq_read_frame(seq, index, &gfit)) {
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

	update_MenuItem();		// initialize menu gui
	display_image_number(index);	// in the gray window
	sequence_list_change_current();
	adjust_exclude(index, FALSE);	// check or uncheck excluded checkbox
	adjust_refimage(index);	// check or uncheck reference image checkbox
	update_used_memory();
	return 0;
}

/**
 * Computes size of an opened sequence in bytes for a passed number of frames.
 * For SER or films, it returns the size of the file.
 * For FITS sequences, the reference image's size is used as baseline.
 * Unsupported for internal sequences.
 * @param seq input sequence
 * @param nb_frames number of frames to compute the size of the sequence of
 * @return the size of the sequence in bytes, or -1 if an error happened.
 */
int64_t seq_compute_size(sequence *seq, int nb_frames) {
	int64_t frame_size, size = -1LL;
	char filename[256];
	GStatBuf sts;
	int ref;

	switch(seq->type) {
	case SEQ_SER:
		size = ser_compute_file_size(seq->ser_file, nb_frames);
		break;
	case SEQ_REGULAR:
		ref = sequence_find_refimage(seq);
		if (fit_sequence_get_image_filename(seq, ref, filename, TRUE)) {
			if (!g_lstat(filename, &sts)) {				
#ifndef _WIN32
				if ( S_ISLNK(sts.st_mode) )
#else
				if (GetFileAttributesA(filename) & FILE_ATTRIBUTE_REPARSE_POINT )
#endif					
				{
					gchar *target_link = g_file_read_link(filename, NULL); 
					if ( !g_lstat(target_link, &sts) )
					{
						frame_size = sts.st_size;       // force 64 bits
					}
					else {
						fprintf(stderr, "Could not open reference image of the sequence\n");
						return size;
					}
					g_free(target_link);
				}
				else
				{
					frame_size = sts.st_size;       // force 64 bits
				}

				size = frame_size * nb_frames;
			}
		}
		break;
#ifdef HAVE_FFMS2
	case SEQ_AVI:
		if (g_stat(seq->film_file->filename, &sts) == 0) {
			// this is a close approximation
			frame_size = sts.st_size / seq->film_file->frame_count;
			size = nb_frames * frame_size;
		}
		break;
#endif
	default:
		fprintf(stderr, "Failure: computing sequence size on internal sequence is unsupported\n");
	}
	return size;
}

/**
 * Check if a sequence with a basename 'basename' already exists
 * @param filename
 * @return TRUE if the name already exists, FALSE otherwise
 */
gboolean check_if_seq_exist(gchar *basename) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		fprintf(stderr, "check_if_seq_exist: %s\n", error->message);
		g_error_free(error);
		g_free(com.wd);
		com.wd = NULL;
		return 1;
	}

	while ((file = g_dir_read_name(dir)) != NULL) {
		gchar *seq = g_strdup_printf("%s.seq", basename);

		if (!g_ascii_strcasecmp(seq, file)) {
			g_free(seq);
			return TRUE;
		}
		g_free(seq);
	}
	return FALSE;
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
			snprintf(name_buf, 255, "%s_%d.ser", seq->seqname,  index);
			name_buf[255] = '\0';
			return name_buf;
#ifdef HAVE_FFMS2
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
#ifdef HAVE_FFMS2
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
	full_stats_invalidation_from_fit(dest);
	copy_seq_stats_to_fit(seq, index, dest);
	return 0;
}

/* same as seq_read_frame above, but creates an image the size of the selection
 * rectangle only. layer is set to the layer number in the read partial frame.
 * The partial image result is only one-channel deep, so it cannot be used to
 * have a partial RGB image. */
int seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry) {
	char filename[256];
#ifdef HAVE_FFMS2
	fits tmp_fit;
#endif
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
			if (ser_read_opened_partial_fits(seq->ser_file, layer, index, dest, area)) {
				siril_log_message(_("Could not load frame %d from SER sequence %s\n"),
						index, seq->seqname);
				return 1;
			}

			break;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			assert(seq->film_file);
			memset(&tmp_fit, 0, sizeof(fits));
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
			PRINT_ALLOC_ERR;
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
					PRINT_ALLOC_ERR;
					return 1;
				}
			}
			if (_allocate_sequence_locks(seq))
				return 1;

			fit_sequence_get_image_filename(seq, index, filename, TRUE);
			siril_fits_open_diskfile(&seq->fptr[index], filename, READONLY, &status);
			if (status) {
				fits_report_error(stderr, status);
				return status;
			}
			/* should we check image parameters here? such as bitpix or naxis */
			break;
		case SEQ_SER:
			assert(seq->ser_file->file == NULL);
			break;
#ifdef HAVE_FFMS2
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
static void set_fwhm_star_as_star_list_with_layer(sequence *seq, int layer) {
	assert(seq->regparam);
	/* we chose here the first layer that has been allocated, which doesn't
	 * mean it contains data for all images. Handle with care. */
	if (seq->regparam && layer >= 0 && layer < seq->nb_layers && seq->regparam[layer] &&
			seq->regparam[layer][seq->current].fwhm_data && !com.stars) {
		com.stars = malloc(2 * sizeof(fitted_PSF *));
		com.stars[0] = seq->regparam[layer][seq->current].fwhm_data;
		com.stars[1] = NULL;
		// this is freed in PSF_list.c:clear_stars_list()
		com.star_is_seqdata = TRUE;
	}
}

// cannot be called in the worker thread
void set_fwhm_star_as_star_list(sequence *seq) {
	int layer = get_registration_layer(seq);
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
	if (seq->fixed <= 1) {
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

char *fit_sequence_get_image_filename_prefixed(sequence *seq, const char *prefix, int index) {
	char format[16];
	GString *str = g_string_sized_new(70);
	sprintf(format, "%%s%%s%%0%dd%%s", seq->fixed);
	g_string_printf(str, format, prefix,
			seq->seqname, seq->imgparam[index].filenum,
			com.ext);
	return g_string_free(str, FALSE);
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

void remove_prefixed_sequence_files(sequence *seq, const char *prefix) {
	int i, len;
	gchar *basename, *seqname;
	if (!prefix || prefix[0] == '\0')
		return;
	basename = g_path_get_basename(seq->seqname);
	len = strlen(basename) + 5 + strlen(prefix);
	seqname = malloc(len);
	g_snprintf(seqname, len, "%s%s.seq", prefix, basename);
	siril_debug_print("Removing %s\n", seqname);
	g_unlink(seqname); // removing the seqfile
	free(seqname);
	g_free(basename);

	switch (seq->type) {
	default:
	case SEQ_REGULAR:
		for (i = 0; i < seq->number; i++) {
			// TODO: use com.cache_upscaled and the current sequence
			// filter to leave the images to be up-scaled.
			char *filename = fit_sequence_get_image_filename_prefixed(
					seq, prefix, i);
			siril_debug_print("Removing %s\n", filename);
			g_unlink(filename);
			free(filename);
		}
		break;
	case SEQ_SER:
		basename = seq->ser_file->filename;
		len = strlen(basename) + strlen(prefix) + 1;
		seqname = malloc(len);
		g_snprintf(seqname, len, "%s%s", prefix, basename);
		siril_debug_print("Removing %s\n", seqname);
		g_unlink(seqname);
		break;
	}
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
	seq->upscale_at_stacking = 1.0;
}

/* call this to close a sequence. Second arg must be FALSE for com.seq
 * WARNING: the data is not reset to NULL, if seq is to be reused,
 * initialize_sequence() must be called on it right after free_sequence()
 * (= do it for com.seq) */
void free_sequence(sequence *seq, gboolean free_seq_too) {
	int layer, j;

	if (seq == NULL) return;
	// free regparam
	if (seq->nb_layers > 0 && seq->regparam) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->regparam[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->regparam[layer][j].fwhm_data &&
							(!seq->photometry[0] ||
							 seq->regparam[layer][j].fwhm_data != seq->photometry[0][j])) // avoid double free
						free(seq->regparam[layer][j].fwhm_data);
				}
				free(seq->regparam[layer]);
			}
		}
		free(seq->regparam);
	}
	// free stats
	if (seq->nb_layers > 0 && seq->stats) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->stats[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->stats[layer][j])
						free_stats(seq->stats[layer][j]);
				}
				free(seq->stats[layer]);
			}
		}
		free(seq->stats);
	}
	// free backup regparam
	if (seq->nb_layers > 0 && seq->regparam_bkp) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->regparam_bkp[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->regparam_bkp[layer][j].fwhm_data &&
							(!seq->photometry[0] ||
							 seq->regparam_bkp[layer][j].fwhm_data != seq->photometry[0][j])) // avoid double free
						free(seq->regparam_bkp[layer][j].fwhm_data);
				}
				free(seq->regparam_bkp[layer]);
			}
		}
		free(seq->regparam_bkp);
	}
	// free backup stats
	if (seq->nb_layers > 0 && seq->stats_bkp) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->stats_bkp[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->stats_bkp[layer][j])
						free_stats(seq->stats_bkp[layer][j]);
				}
				free(seq->stats_bkp[layer]);
			}
		}
		free(seq->stats_bkp);
	}

	// free name of the layers
	if (seq->nb_layers > 0) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			free(seq->layers[layer].name);
		}
	}

	for (j=0; j<seq->number; j++) {
		if (seq->fptr && seq->fptr[j]) {
			int status = 0;
			fits_close_file(seq->fptr[j], &status);
		}
		if (seq->imgparam) {
			if (seq->imgparam[j].date_obs)
				free(seq->imgparam[j].date_obs);
		}
	}
	if (seq->seqname)	free(seq->seqname);
	if (seq->layers)	free(seq->layers);
	if (seq->imgparam)	free(seq->imgparam);
	if (seq->fptr)		free(seq->fptr);

#ifdef _OPENMP
	if (seq->fd_lock) {
		for (j=0; j<seq->number; j++) {
			omp_destroy_lock(&seq->fd_lock[j]);
		}
		free(seq->fd_lock);
	}
#endif

	if (seq->ser_file) {
		ser_close_file(seq->ser_file);	// frees the data too
		free(seq->ser_file);
	}
#ifdef HAVE_FFMS2
	if (seq->film_file) {
		film_close_file(seq->film_file);	// frees the data too
		free(seq->film_file);
	}
#endif

	if (seq->internal_fits) {
		// Compositing still uses references to the images in the sequence
		//for (j=0; j<seq->number; j++)
		//	clearfits(seq->internal_fits[j]);
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

	for (j = 0; j < MAX_SEQPSF && seq->photometry[j]; j++) {
		free_photometry_set(seq, j);
	}

	if (free_seq_too)	free(seq);
}

gboolean sequence_is_loaded() {
	return (com.seq.seqname != NULL && com.seq.imgparam != NULL);
}

/* Close the com.seq sequence */
void close_sequence(int loading_another) {
	fprintf(stdout, "MODE: closing sequence\n");
	if (sequence_is_loaded()) {
		siril_log_message(_("Closing sequence %s\n"), com.seq.seqname);
		if (!com.headless) {
			free_cbbt_layers();
			clear_sequence_list();
		}
		if (com.seq.needs_saving)
			writeseqfile(&com.seq);
		free_sequence(&com.seq, FALSE);
		initialize_sequence(&com.seq, FALSE);
		if (!com.headless) {
			clear_stars_list();
			update_stack_interface(TRUE);
		}
		if (!loading_another && !com.headless) {
			// unselect the sequence in the sequence list
			GtkComboBox *seqcombo = GTK_COMBO_BOX(lookup_widget("sequence_list_combobox"));
			gtk_combo_box_set_active(seqcombo, -1);
		}
	}
}

/* if no reference image has been set, return the index of an image that is
 * selected in the sequence, the best of the first registration data found if
 * any, the first otherwise */
int sequence_find_refimage(sequence *seq) {
	if (seq->reference_image != -1)
		return seq->reference_image;
	if (seq->type == SEQ_INTERNAL)
		return 1; // green channel
	int layer, image, best = -1;
	for (layer = 0; layer < seq->nb_layers; layer++) {
		if (seq->regparam[layer]) {
			gboolean use_fwhm;
			double best_val;
			if (seq->regparam[layer][0].fwhm > 0.0) {
				use_fwhm = TRUE;
				best_val = 1000000.0;
			} else if (seq->regparam[layer][0].quality > 0.0) {
				use_fwhm = FALSE;
				best_val = 0.0;
			}
			else continue;

			for (image = 0; image < seq->number; image++) {
				if (!seq->imgparam[image].incl)
					continue;
				if (use_fwhm && seq->regparam[layer][image].fwhm > 0 &&
						seq->regparam[layer][image].fwhm < best_val) {
					best_val = seq->regparam[layer][image].fwhm;
					best = image;
				} else if (seq->regparam[layer][image].quality > 0 &&
						seq->regparam[layer][image].quality > best_val) {
					best_val = seq->regparam[layer][image].quality;
					best = image;
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

/* assign shift values for registration data of a sequence, depending on its type and sign */
void set_shifts(sequence *seq, int frame, int layer, float shiftx, float shifty, gboolean data_is_top_down) {
	if (seq->regparam[layer]) {
		seq->regparam[layer][frame].shiftx = shiftx;
		seq->regparam[layer][frame].shifty = data_is_top_down ? -shifty : shifty;
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

fits *internal_sequence_get(sequence *seq, int index) {
	if (index > seq->number)
		return NULL;
	return seq->internal_fits[index];
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

	stop_processing_thread();
	if (!args->retvalue) {
		char *rseqname = malloc(
				strlen(args->prefix) + strlen(args->seq->seqname) + 5);

		sprintf(rseqname, "%s%s.seq", args->prefix, args->seq->seqname);
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
			PRINT_ALLOC_ERR;
			args->retvalue = 1;
			siril_add_idle(end_crop_sequence, args);
		}
		sprintf(dest, "%s%s.ser", args->prefix, args->seq->seqname);
		if (ser_create_file(dest, ser_file, TRUE, com.seq.ser_file)) {
			siril_log_message(_("Creating the SER file failed, aborting.\n"));
			free(ser_file);
			ser_file = NULL;
			args->retvalue = 1;
			siril_add_idle(end_crop_sequence, args);
		}
	}
	set_progress_bar_data(_("Processing..."), PROGRESS_RESET);
	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		if (!get_thread_run() || args->retvalue)
			break;
		fits fit;
		memset(&fit, 0, sizeof(fits));
		ret = seq_read_frame(args->seq, frame, &fit);
		if (!ret) {
			char dest[256], filename[256];

			crop(&fit, &(args->area));
			switch (args->seq->type) {
			case SEQ_REGULAR:
				fit_sequence_get_image_filename(args->seq, frame, filename,
				TRUE);
				sprintf(dest, "%s%s", args->prefix, filename);
				args->retvalue = savefits(dest, &fit);
				break;
			case SEQ_SER:
				if (ser_file) {
					ser_file->image_width = fit.rx;
					ser_file->image_height = fit.ry;
					if (ser_write_frame_from_fit(ser_file, &fit, frame)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						args->retvalue = 1;
					}
				}
				break;
			default:
				args->retvalue = 1;	// should not happen
			}

			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / args->seq->number);
		}
		clearfits(&fit);
	}
	if (args->seq->type == SEQ_SER) {
		ser_write_and_close(ser_file);
		free(ser_file);
	}
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	siril_add_idle(end_crop_sequence, args);
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

	/* The seqpsf result for each image, list of seqpsf_data */
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
int seqpsf_image_hook(struct generic_seq_args *args, int out_index, int index, fits *fit, rectangle *area) {
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;
	struct seqpsf_data *data = malloc(sizeof(struct seqpsf_data));
	data->image_index = index;

	rectangle psfarea = { .x = 0, .y = 0, .w = fit->rx, .h = fit->ry };
	data->psf = psf_get_minimisation(fit, 0, &psfarea, !spsfargs->for_registration, TRUE);
	if (data->psf) {
		data->psf->xpos = data->psf->x0 + area->x;
		if (fit->top_down)
			data->psf->ypos = data->psf->y0 + area->y;
		else data->psf->ypos = area->y + area->h - data->psf->y0;

		/* let's move args->area to center it on the star */
		if (spsfargs->framing == FOLLOW_STAR_FRAME) {
			args->area.x = round_to_int(data->psf->xpos - args->area.w*0.5);
			args->area.y = round_to_int(data->psf->ypos - args->area.h*0.5);
			//fprintf(stdout, "moving area to %d, %d\n", args->area.x, args->area.y);
		}

		if (!args->seq->imgparam[index].date_obs && fit->date_obs[0] != '\0')
			args->seq->imgparam[index].date_obs = strdup(fit->date_obs);
		data->exposure = fit->exposure;
	}
	else {
		if (spsfargs->framing == FOLLOW_STAR_FRAME) {
			siril_log_color_message(_("No star found in the area image %d around %d,%d"
						" (use a larger area?)\n"),
					"red", index, area->x, area->y);
		} else {
			siril_log_color_message(_("No star found in the area image %d around %d,%d"
					" (use 'follow star' option?)\n"),
				"red", index, area->x, area->y);
		}
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
	int photometry_index = 0;
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

	GSList *iterator;
	for (iterator = spsfargs->list; iterator; iterator = iterator->next) {
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
			if (data->psf) {
				seq->regparam[layer][data->image_index].fwhm = data->psf->fwhmx;
				seq->regparam[layer][data->image_index].roundness =
					data->psf->fwhmy / data->psf->fwhmx;
			}
		}

		// for photometry use: store data in seq->photometry
		if (!spsfargs->for_registration) {
			seq->photometry[photometry_index][data->image_index] = data->psf;
		}
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
			fill_sequence_list(seq, layer, FALSE);

		set_layers_for_registration();	// update display of available reg data
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
	}

proper_ending:
	dont_stop_thread = args->already_in_a_thread;
	if (spsfargs->list)
		g_slist_free(spsfargs->list);
	free(spsfargs);
	adjust_sellabel();

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
int seqpsf(sequence *seq, int layer, gboolean for_registration, gboolean regall,
		framing_mode framing, gboolean run_in_thread) {

	if (framing == REGISTERED_FRAME && !seq->regparam[layer])
		framing = ORIGINAL_FRAME;

	if (com.selection.w <= 0 || com.selection.h <= 0){
		siril_log_message(_("Select an area first\n"));
		return 1;
	}
	if (framing == FOLLOW_STAR_FRAME)
		siril_log_color_message(_("The sequence analysis of the PSF will use a sliding selection area centred on the previous found star; this disables parallel processing.\n"), "salmon");
	else if (framing == REGISTERED_FRAME)
		siril_log_color_message(_("The sequence analysis of the PSF will use registration data to move the selection area for each image; this is compatible with parallel processing.\n"), "salmon");

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
	args->filtering_criterion = (regall == TRUE) ? seq_filter_all : seq_filter_included;
	args->nb_filtered_images = (regall == TRUE) ? seq->number : seq->selnum;
	args->prepare_hook = NULL;
	args->finalize_hook = NULL;
	args->image_hook = seqpsf_image_hook;
	args->idle_function = end_seqpsf;
	args->stop_on_error = FALSE;
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

