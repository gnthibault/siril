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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "core/siril.h"
#include "algos/statistics.h"
#include "io/ser.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#ifdef HAVE_FFMS2
#include "io/films.h"
#endif

/* seqfile version history *
 * no version up to 0.9.9
 * version 1 introduced roundness in regdata, 0.9.9
 * version 2 allowed regdata to be stored for CFA SER sequences, 0.9.11
 */
#define CURRENT_SEQFILE_VERSION 2	// to increment on format change

/* File format (lines starting with # are comments, lines that are (for all
 * something) need to be in all in sequence of this only type of line):
 *
 * S sequence_name beg number selnum fixed reference_image [version]
 * L nb_layers
 * (for all images) I filenum incl [stats+] <- stats added at some point, removed in 0.9.9
 * (for all layers (x)) Rx regparam+
 * TS | TA (type for ser or film)
 * U up-scale_ratio
 * (for all images (y) and layers (x)) Mx-y stats+
 */

/* name is sequence filename, with or without .seq extension
 * It should always be used with seq_check_basic_data() because on first loading
 * of a .seq that was created from scan of the filesystem, number of layers and
 * image size are unknown and some properties of the sequence are null or unset.
 * Returns NULL if the sequence could not be loaded.
 */
sequence * readseqfile(const char *name){
	char line[512], *scanformat;
	char filename[512], *seqfilename;
	int i, nbsel, nb_tokens, allocated = 0, current_layer = -1, image;
	int to_backup = 0, version = -1;
	FILE *seqfile;
       	sequence *seq;
	imstats *stats;
	regdata *regparam;

	if (!name) return NULL;
	fprintf(stdout, "Reading sequence file `%s'.\n", name);

	if(!ends_with(name, ".seq")){
		seqfilename = malloc(strlen(name) + 6);	/* 6 stands for a max length of 4 + '.' + '\0' */
		sprintf(seqfilename, "%s.seq", name);
	} else {
		seqfilename = strdup(name);
	}

	if ((seqfile = g_fopen(seqfilename, "r")) == NULL) {
		fprintf(stderr, "Reading sequence failed, file cannot be opened: %s.\n", seqfilename);
		free(seqfilename);
		return NULL;
	}

	seq = calloc(1, sizeof(sequence));
	initialize_sequence(seq, TRUE);
	i=0;
	while (fgets(line, 511, seqfile)) {
		switch (line[0]) {
			case '#':
				continue;
			case 'S':
				/* The double quote as sequence name is a sequence with no name.
				 * Such sequences don't exist anymore. */
				assert(line[2] != '"');
				if (line[2] == '\'')	/* new format, quoted string */
					scanformat = "'%511[^']' %d %d %d %d %d %d";
				else scanformat = "%511s %d %d %d %d %d %d";

				if(sscanf(line+2, scanformat,
							filename, &seq->beg, &seq->number,
							&seq->selnum, &seq->fixed,
							&seq->reference_image, &version) < 6 ||
						allocated != 0){
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				if (seq->number == 0) {
					fprintf(stderr, "readseqfile: sequence is empty?\n");
					goto error;
				}
				if (version > CURRENT_SEQFILE_VERSION)
					siril_log_message(_("This sequence file was created by a version of "
								"siril that is newer than this one, it may not "
								"be loaded as expected\n"),
							"salmon");
				/* for now, only the R* line is not supported in the previous version */
				seq->seqname = strdup(filename);
				seq->imgparam = calloc(seq->number, sizeof(imgdata));
				allocated = 1;
				break;

			case 'L':
				/* for now, the L line stores the number of layers for each image. */
				if (line[1] == ' ') {
					int nbl_backup = seq->nb_layers;
					if (sscanf(line+2, "%d", &seq->nb_layers) != 1) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
					/* seq->nb_layers can be -1 when the sequence has not been
					 * opened for the first time, but it may already have been
					 * set in SER opening below, so we keep the backup in this
					 * case */
					if (nbl_backup > 0 && ser_is_cfa(seq->ser_file)) {
						if (com.debayer.open_debayer)
							seq->nb_layers = nbl_backup;
						else seq->nb_layers = 1;
					}
					// else if nbl_backup is 3 but opening debayer is not
					// enabled, we keep 1 in the nb_layers, which will be set in
					// the seq_check_basic_data() call later
					if (seq->nb_layers >= 1) {
						seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
						seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
						if (ser_is_cfa(seq->ser_file))
							seq->regparam_bkp = calloc(3, sizeof(regdata*));
					}
				} else if (line[1] >= '0' && line[1] <= '9') {
					/* in the future, wavelength and name of each layer will be added here */
				}
				break;

			case 'I':
				/* First sequence file format was I filenum and incl.
				 * A later sequence file format added the stats to this.
				 * The current file format comes back to the first,
				 * moving stats to the M-line. */
				stats = NULL;
				if (!seq->imgparam) {
					fprintf(stderr, "readseqfile: sequence file format error, missing S line\n");
					goto error;
				}
				allocate_stats(&stats);
				nb_tokens = sscanf(line + 2,
						"%d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
						&(seq->imgparam[i].filenum),
						&(seq->imgparam[i].incl),
						&(stats->mean),
						&(stats->median),
						&(stats->sigma),
						&(stats->avgDev),
						&(stats->mad),
						&(stats->sqrtbwmv),
						&(stats->location),
						&(stats->scale),
						&(stats->min),
						&(stats->max));
				if (nb_tokens == 12) {
					add_stats_to_seq(seq, i, 0, stats);
					free_stats(stats);	// we unreference it here
				} else {
					free_stats(stats);
					if (nb_tokens != 2) {
						fprintf(stderr, "readseqfile: sequence file format error: %s\n", line);
						goto error;
					}
				}
				++i;
				break;

			case 'R':
				/* registration info */
				if (line[1] == '*') {
					/* these are registration data for the CFA channel, the
					 * star is a way to differentiate stats belonging to
					 * CFA and those belonging to the demosaiced red
					 * channel, both would have layer number 0 otherwise */
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.debayer.open_debayer) {
						siril_debug_print("- using CFA registration info\n");
						to_backup = 0;
					} else {
						siril_debug_print("- backing up CFA registration info\n");
						to_backup = 1;
					}
					current_layer = 0;
				}
				else {
					to_backup = 0;
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.debayer.open_debayer) {
						to_backup = 1;
						siril_debug_print("- stats: backing up demosaiced stats\n");
					}
					current_layer = line[1] - '0';
				}

				if (current_layer < 0 || current_layer > 9 ||
						(!seq->cfa_opened_monochrome && current_layer >= seq->nb_layers) ||
						(seq->cfa_opened_monochrome && current_layer >= 3)) {
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}

				if (to_backup)
					regparam = seq->regparam_bkp[current_layer];
				else regparam = seq->regparam[current_layer];

				if (!regparam) {
					regparam = calloc(seq->number, sizeof(regdata));
					if (!regparam) {
						PRINT_ALLOC_ERR;
						goto error;
					}
					i = 0;	// one line per image, starting with 0
					// reassign, because we didn't use a pointer
					if (to_backup)
						seq->regparam_bkp[current_layer] = regparam;
					else seq->regparam[current_layer] = regparam;
				}
				if (i >= seq->number) {
					fprintf(stderr, "\nreadseqfile ERROR: out of array bounds in reg info!\n\n");
					goto error;
				}
				if (version < 1) {
					float rot_centre_x, rot_centre_y, angle;
					nb_tokens = sscanf(line+3, "%f %f %g %g %g %g %lg",
							&(regparam[i].shiftx),
							&(regparam[i].shifty),
							&rot_centre_x, &rot_centre_y,
							&angle,
							&(regparam[i].fwhm),
							&(regparam[i].quality));
					if (nb_tokens != 7) {
						if (nb_tokens == 3) {
							// old format, with quality as third token
							regparam[i].quality = rot_centre_x;
							// the rest is already zero due to the calloc
						} else {
							fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
							goto error;
						}
					}
				} else {
					// new file format with roundness instead of weird things
					if (sscanf(line+3, "%f %f %g %g %lg",
								&(regparam[i].shiftx),
								&(regparam[i].shifty),
								&(regparam[i].fwhm),
								&(regparam[i].roundness),
								&(regparam[i].quality)) != 5) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
				}
				++i;
				break;

			case 'T':
				/* sequence type (several files or a single file) */
				if (line[1] == 'S') {
					seq->type = SEQ_SER;
#ifdef HAVE_FFMS2
					seq->ext = "ser";
#endif
					if (seq->ser_file) break;
					seq->ser_file = malloc(sizeof(struct ser_struct));
					ser_init_struct(seq->ser_file);
					seqfilename[strlen(seqfilename)-1] = 'r';
					if (ser_open_file(seqfilename, seq->ser_file)) {
						free(seq->ser_file);
						seq->ser_file = NULL;
						goto error;
					}
					ser_display_info(seq->ser_file);

					if (ser_is_cfa(seq->ser_file)) {
						if (!com.debayer.open_debayer) {
							// we set this flag instead of relying on the
							// com.debayer.open_debayer flag which varies
							// as the user changes the GUI
							seq->cfa_opened_monochrome = TRUE;
						}
						seq->nb_layers = 3;
						if (seq->regparam)
							seq->regparam = realloc(seq->regparam, seq->nb_layers * sizeof(regdata *));
						if (seq->layers)
							seq->layers = realloc(seq->regparam, seq->nb_layers * sizeof(layer_info));
						seq->needs_saving = TRUE;
					}
				}
#ifdef HAVE_FFMS2
				else if (line[1] == 'A') {
					seq->type = SEQ_AVI;
					if (seq->film_file) break;
					seq->film_file = malloc(sizeof(struct film_struct));
					int ii = 0, nb_film = get_nb_film_ext_supported();

					gchar *filmname = NULL;
					while (ii < nb_film) {
						GString *filmString;
						/* test for extension in lowercase */
						filmString = g_string_new(seqfilename);
						filmString = g_string_truncate(filmString, strlen(seqfilename) - 3);
						filmString = g_string_append(filmString, supported_film[i].extension);
						g_free(filmname);
						filmname = g_string_free(filmString, FALSE);

						if (g_file_test(filmname, G_FILE_TEST_EXISTS)) {
							break;
						} else {
							int len_ext;
							gchar *upcase;

							g_free(filmname);
							filmname = NULL;

							filmString = g_string_new(seqfilename);
							g_string_truncate(filmString, strlen(seqfilename) - 3);
							len_ext = strlen(supported_film[i].extension);
							upcase = g_ascii_strup(supported_film[i].extension, len_ext);
							filmString = g_string_append(filmString, upcase);
							filmname = g_string_free(filmString, FALSE);
							g_free(upcase);

							if (g_file_test(filmname, G_FILE_TEST_EXISTS)) {
								break;
							}
						}
						ii++;
					}

					if (film_open_file(filmname, seq->film_file)) {
						free(seq->film_file);
						seq->film_file = NULL;
						g_free(filmname);
						goto error;
					}
					else {
						film_display_info(seq->film_file);
						seq->ext = strdup(get_filename_ext(seq->film_file->filename));
						g_free(filmname);
					}
				}
				else seq->ext = "fit";
#endif
				break;

			case 'U':
				/* up-scale factor for stacking. Used in simplified stacking for
				 * shift-only registrated sequences, up-scale will be done at
				 * stack-time. */
				if (line[1] == ' ' &&
						sscanf(line+2, "%lg", &seq->upscale_at_stacking) != 1) {
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				break;
			case 'M':
				/* stats may not exist for all images and layers so we use
				 * indices for them, the line is Mx-y with x the layer number
				 * and y the image index */

				if (line[1] == '*') {
					/* these are stats for the CFA channel, the star is a
					 * way to differentiate stats belonging to CFA and
					 * those belonging to the demosaiced red channel, both
					 * would have layer number 0 otherwise */
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.debayer.open_debayer) {
						siril_debug_print("- stats: using CFA stats\n");
						to_backup = 0;
					} else {
						siril_debug_print("- stats: backing up CFA stats\n");
						to_backup = 1;
					}
					current_layer = 0;
				}
				else {
					to_backup = 0;
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.debayer.open_debayer) {
						to_backup = 1;
						siril_debug_print("- stats: backing up demosaiced stats\n");
					}
					current_layer = line[1] - '0';
				}

				if (current_layer < 0 || current_layer > 9 || line[2] != '-') {
					fprintf(stderr, "readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				/*if (current_layer >= seq->nb_layers) {
				// it may happen when opening a CFA file in monochrome
				break;
				}*/
				stats = NULL;
				allocate_stats(&stats);
				nb_tokens = sscanf(line + 3,
						"%d %ld %ld %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
						&image,
						&(stats->total),
						&(stats->ngoodpix),
						&(stats->mean),
						&(stats->median),
						&(stats->sigma),
						&(stats->avgDev),
						&(stats->mad),
						&(stats->sqrtbwmv),
						&(stats->location),
						&(stats->scale),
						&(stats->min),
						&(stats->max),
						&(stats->normValue),
						&(stats->bgnoise));
				if (nb_tokens == 15) {
					if (to_backup)
						add_stats_to_seq_backup(seq, image, current_layer, stats);
					else add_stats_to_seq(seq, image, current_layer, stats);
					free_stats(stats);	// we unreference it here
				} else {
					free_stats(stats);
					fprintf(stderr, "readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				break;
		}
	}
	if (!allocated) {
		siril_log_message(_("The sequence file %s seems to be corrupted\n"), seqfilename);
		goto error;
	}
	seq->needs_saving = FALSE;	// loading stats sets it to true
	fclose(seqfile);
	seq->end = seq->imgparam[seq->number-1].filenum;
	seq->current = -1;
	for (i=0, nbsel=0; i<seq->number; i++)
		if (seq->imgparam[i].incl)
			nbsel++;
	if (nbsel != seq->selnum) {
		siril_log_message(_("Fixing the selection number in the .seq file (%d) to the actual value (%d) (not saved)\n"), seq->selnum, nbsel);
		seq->selnum = nbsel;
	}
	
	// copy some regparam_bkp to regparam if it applies
	if (ser_is_cfa(seq->ser_file) && com.debayer.open_debayer &&
			seq->regparam_bkp && seq->regparam_bkp[0] &&
			seq->regparam && seq->nb_layers == 3 && !seq->regparam[1]) {
		siril_log_color_message(_("%s: Copying registration data from non-demosaiced layer to green layer\n"), "salmon", seqfilename);
		seq->regparam[1] = calloc(seq->number, sizeof(regdata));
		for (image = 0; image < seq->number; image++) {
			memcpy(&seq->regparam[1][image], &seq->regparam_bkp[0][image], sizeof(regdata));
		}
	}

	update_used_memory();
	free(seqfilename);
	return seq;
error:
	fclose(seqfile);
	if (seq->seqname)
		free(seq->seqname);
	free(seq);
	siril_log_message(_("Could not load sequence %s\n"), name);
	update_used_memory();
	free(seqfilename);
	return NULL;
}

/* Saves the sequence in the seqname.seq file. */
int writeseqfile(sequence *seq){
	char *filename;
	FILE *seqfile;
	int i, layer;

	if (!seq->seqname || seq->seqname[0] == '\0') return 1;
	filename = malloc(strlen(seq->seqname)+5);
	sprintf(filename, "%s.seq", seq->seqname);
	seqfile = g_fopen(filename, "w+t");
	if (seqfile == NULL) {
		perror("writeseqfile, fopen");
		fprintf(stderr, "Writing sequence file: cannot open %s for writing\n", filename);
		free(filename);
		return 1;
	}
	fprintf(stdout, "Writing sequence file %s\n", filename);
	free(filename);

	fprintf(seqfile,"#Siril sequence file. Contains list of files (images), selection, and registration data\n");
	fprintf(seqfile,"#S 'sequence_name' start_index nb_images nb_selected fixed_len reference_image version\n");
	fprintf(stderr,"S '%s' %d %d %d %d %d %d\n", 
			seq->seqname, seq->beg, seq->number, seq->selnum, seq->fixed, seq->reference_image, CURRENT_SEQFILE_VERSION);
	fprintf(seqfile,"S '%s' %d %d %d %d %d %d\n", 
			seq->seqname, seq->beg, seq->number, seq->selnum, seq->fixed, seq->reference_image, CURRENT_SEQFILE_VERSION);
	if (seq->type != SEQ_REGULAR) {
		/* sequence type, not needed for regular, S for ser, A for avi */
		fprintf(stderr, "T%c\n", seq->type == SEQ_SER ? 'S' : 'A');
		fprintf(seqfile, "T%c\n", seq->type == SEQ_SER ? 'S' : 'A');
	}

	if (seq->upscale_at_stacking != 1.0) {
		// until we have a real drizzle
		fprintf(stderr, "U %g\n", seq->upscale_at_stacking);
		fprintf(seqfile, "U %g\n", seq->upscale_at_stacking);
	}

	fprintf(stderr, "L %d\n", seq->nb_layers);
	fprintf(seqfile, "L %d\n", seq->nb_layers);

	for(i=0; i < seq->number; ++i){
		fprintf(seqfile,"I %d %d\n",
				seq->imgparam[i].filenum, 
				seq->imgparam[i].incl);
	}

	for (layer = 0; layer < seq->nb_layers; layer++) {
		if (seq->regparam[layer]) {
			for (i=0; i < seq->number; ++i) {
				fprintf(seqfile, "R%c %f %f %g %g %g\n",
						seq->cfa_opened_monochrome ? '*' : '0' + layer,
						seq->regparam[layer][i].shiftx,
						seq->regparam[layer][i].shifty,
						seq->regparam[layer][i].fwhm,
						seq->regparam[layer][i].roundness,
						seq->regparam[layer][i].quality
				       );
			}
		}
		if (seq->stats && seq->stats[layer]) {
			for (i=0; i < seq->number; ++i) {
				if (!seq->stats[layer][i]) continue;

				fprintf(seqfile, "M%c-%d %ld %ld %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
						seq->cfa_opened_monochrome ? '*' : '0' + layer, i,
						seq->stats[layer][i]->total,
						seq->stats[layer][i]->ngoodpix,
						seq->stats[layer][i]->mean,
						seq->stats[layer][i]->median,
						seq->stats[layer][i]->sigma,
						seq->stats[layer][i]->avgDev,
						seq->stats[layer][i]->mad,
						seq->stats[layer][i]->sqrtbwmv,
						seq->stats[layer][i]->location,
						seq->stats[layer][i]->scale,
						seq->stats[layer][i]->min,
						seq->stats[layer][i]->max,
						seq->stats[layer][i]->normValue,
						seq->stats[layer][i]->bgnoise);

			}
		}
	}
	for (layer = 0; layer < 3; layer++) {
		if (seq->regparam_bkp && seq->regparam_bkp[layer]) {
			for (i=0; i < seq->number; ++i) {
				fprintf(seqfile, "R%c %f %f %g %g %g\n",
						seq->cfa_opened_monochrome ? '0' + layer : '*',
						seq->regparam_bkp[layer][i].shiftx,
						seq->regparam_bkp[layer][i].shifty,
						seq->regparam_bkp[layer][i].fwhm,
						seq->regparam_bkp[layer][i].roundness,
						seq->regparam_bkp[layer][i].quality
				       );
			}
		}
		if (seq->stats_bkp && seq->stats_bkp[layer]) {
			for (i=0; i < seq->number; ++i) {
				if (!seq->stats_bkp[layer][i]) continue;

				fprintf(seqfile, "M%c-%d %ld %ld %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
						seq->cfa_opened_monochrome ? '0' + layer : '*', i,
						seq->stats_bkp[layer][i]->total,
						seq->stats_bkp[layer][i]->ngoodpix,
						seq->stats_bkp[layer][i]->mean,
						seq->stats_bkp[layer][i]->median,
						seq->stats_bkp[layer][i]->sigma,
						seq->stats_bkp[layer][i]->avgDev,
						seq->stats_bkp[layer][i]->mad,
						seq->stats_bkp[layer][i]->sqrtbwmv,
						seq->stats_bkp[layer][i]->location,
						seq->stats_bkp[layer][i]->scale,
						seq->stats_bkp[layer][i]->min,
						seq->stats_bkp[layer][i]->max,
						seq->stats_bkp[layer][i]->normValue,
						seq->stats_bkp[layer][i]->bgnoise);

			}
		}
	}

	fclose(seqfile);
	seq->needs_saving = FALSE;
	return 0;
}

gboolean existseq(const char *name){
	char *filename;
	GStatBuf sts;
	if (!name || name[0] == '\0') return FALSE;
	filename = malloc(strlen(name)+5);
	sprintf(filename, "%s.seq", name);
	if(g_stat(filename, &sts)==0){
		free(filename);
		return TRUE;
	}
	free(filename);
	return FALSE;
}

/* try to create the sequence file for the newly found sequence */
int buildseqfile(sequence *seq, int force_recompute) {
	image_type imagetype;
	int i;
	char *filename;
	imgdata *oldparam;

	if (seq->end <= 0 || !seq->seqname || seq->seqname[0] == '\0') return 1;
	if (existseq(seq->seqname) && !force_recompute) {
		fprintf(stderr,"seqfile '%s.seq' already exists, not recomputing\n", seq->seqname);
		return 0;
	}

	if (force_recompute) {
		for (i = 0; i < seq->nb_layers; i++)
			clear_stats(seq, i);
	}

	filename = malloc(strlen(seq->seqname) + 20);
	if (filename == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	if (seq->type == SEQ_REGULAR) {
		get_possible_image_filename(seq, seq->beg, filename);
		// check if the sequence begins at first_index

		if (stat_file(filename, &imagetype, NULL) || imagetype != TYPEFITS) {
			siril_log_message(_("The sequence %s doesn't start at the frame number %d"
					" with the specified fixed size index (%d). Cannot load.\n"),
					seq->seqname, seq->beg, seq->fixed);
			free(filename);
			return 1;
		}
	}

	int alloc_size = 30;
	//seq->number = 0;
	// fill in one pass: realloc needed
	if (seq->type == SEQ_REGULAR) {
		if (seq->end - seq->beg < 111)
			alloc_size = seq->end - seq->beg + 1;	// last index IS included
	} else alloc_size = seq->end - seq->beg + 1;		// always continuous
	oldparam = seq->imgparam;
	if ((seq->imgparam = realloc(seq->imgparam, alloc_size*sizeof(imgdata))) == NULL) {
		fprintf(stderr, "Could not reallocate image parameters structure in sequence\n");
		if (oldparam) free(oldparam);
		free(filename);
		return 2;
	}
	for (i = seq->beg; i <= seq->end; i++) {
		if (seq->type == SEQ_REGULAR) {
			get_possible_image_filename(seq, i, filename);
			if (!stat_file(filename, &imagetype, NULL) && imagetype == TYPEFITS) {
				if (seq->number+1 > alloc_size-1) {
					alloc_size += 25;
					oldparam = seq->imgparam;
					if (!(seq->imgparam = realloc(seq->imgparam, alloc_size*sizeof(imgdata)))) {
						PRINT_ALLOC_ERR;
						if (oldparam) free(oldparam);
						free(filename);
						return 2;
					}
				}
				seq->imgparam[seq->number].filenum = i;
				seq->imgparam[seq->number].incl = SEQUENCE_DEFAULT_INCLUDE;
				seq->imgparam[seq->number].date_obs = NULL;
				seq->number++;
			}
		} else {
			seq->imgparam[i].filenum = i;
			seq->imgparam[i].incl = SEQUENCE_DEFAULT_INCLUDE;
			seq->imgparam[i].date_obs = NULL;
		}
	}
#if SEQUENCE_DEFAULT_INCLUDE == TRUE
	seq->selnum = seq->number;
#else
	seq->selnum = 0;
#endif
	writeseqfile(seq);

	fprintf(stdout, "Sequence found: %s %d->%d\n", seq->seqname, seq->beg, seq->end);
	free(filename);
	return 0;
}

