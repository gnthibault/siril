/*
 *  match: a package to match lists of stars (or other items)
 *  Copyright (C) 2000  Michael William Richmond
 *
 *  Contact: Michael William Richmond
 *           Physics Department
 *           Rochester Institute of Technology
 *           85 Lomb Memorial Drive
 *           Rochester, NY  14623-5603
 *           E-mail: mwrsps@rit.edu
 *
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/siril_world_cs.h"
#include "misc.h"
#include "project_coords.h"
#include "degtorad.h"
#define RADtoASEC (3600.*180./M_PI)

enum COORD_COL {
	RACOL = 1, DECCOL = 2
};

#undef DEBUG           /* get some of diagnostic output */

static int
proc_star_file(const char *file, int racol, int deccol, double ra, double dec,
		GFile *out_fp, int doASEC);

int convert_catalog_coords(const char *fileA, SirilWorldCS *world_cs, GFile *out) {
	int doASEC = 1;
	double ra = siril_world_cs_get_alpha(world_cs);
	double dec = siril_world_cs_get_delta(world_cs);

	/* now walk through the file and do the dirty work */
	if (proc_star_file(fileA, RACOL, DECCOL, ra, dec, out, doASEC) != SH_SUCCESS) {
		shError("can't process data from file %s", fileA);
		exit(1);
	}

	return (0);
}

/****************************************************************************
 * ROUTINE: proc_star_file
 *
 * walk through the given file, one line at a time.  
 *
 * If the line starts with COMMENT_CHAR, place it into the output stream.
 * If the line is completely blank, place it into the output stream.
 *
 * Otherwise, 
 *   - read in the entire line, 
 *   - figure out the RA and Dec coords of the star
 *   - calculate the (xi, eta) coords of the star projected onto the
 *         plane tangent to the sky at (central_ra, central_dec)
 *   - if the "doASEC" arg is non-zero, convert the values of (xi, eta)
 *         from radians to arcseconds
 *   - print out the line, replacing the "RA" with xi, the "Dec" with eta
 *
 * RETURNS:
 *   SH_SUCCESS            if all goes well
 *   SH_GENERIC_ERROR      if not
 */

static int proc_star_file(const char *file, /* I: name of input file with star list */
int racol, /* I: position of column with RA positions */
int deccol, /* I: position of column with Dec positions */
double central_ra, /* I: central RA of tangent plane (degrees) */
double central_dec, /* I: central Dec of tangent plane (degrees) */
GFile *file_out, /* I: place output into this stream */
int doASEC /* I: if > 0, write offsets in arcsec */
) {
	gchar *line;
	GError *error = NULL;
	GInputStream *input_stream = NULL;
	char col[MAX_DATA_COL + 1][MAX_COL_LENGTH + 1];
	int i, ncol;
	int last_column = -1;
	double raval, decval;
	double /*cent_ra_rad, */cent_dec_rad;
	double dec_rad;
	double delta_ra;
	double xx, yy, xi, eta;

	GFile *file_in = g_file_new_for_path(file);

	input_stream = (GInputStream*) g_file_read(file_in, NULL, &error);

	if (input_stream == NULL) {
		if (error != NULL) {
			shError("proc_star_file: can't open file %s for input", file);
			g_clear_error(&error);
		}

		g_object_unref(file_in);
		return SH_GENERIC_ERROR;
	}

	last_column = (racol > deccol ? racol : deccol);
	//cent_ra_rad = central_ra*DEGTORAD;
	cent_dec_rad = central_dec * DEGTORAD;

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	gsize length = 0;
	while ((line = g_data_input_stream_read_line_utf8(data_input, &length, NULL, NULL))) {

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}

		shAssert(MAX_DATA_COL >= 20);
		ncol = sscanf(line, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
						&(col[0][0]), &(col[1][0]), &(col[2][0]), &(col[3][0]),
						&(col[4][0]), &(col[5][0]), &(col[6][0]), &(col[7][0]),
						&(col[8][0]), &(col[9][0]), &(col[10][0]),
						&(col[11][0]), &(col[12][0]), &(col[13][0]),
						&(col[14][0]), &(col[15][0]), &(col[16][0]),
						&(col[17][0]), &(col[18][0]), &(col[19][0]),
						&(col[20][0]), &(col[21][0]), &(col[22][0]),
						&(col[23][0]), &(col[24][0]), &(col[25][0]),
						&(col[26][0]), &(col[27][0]), &(col[28][0]),
						&(col[29][0]));
		if (last_column > ncol) {
			shError(
					"proc_star_file: not enough entries in following line; skipping");
			shError("  %s", line);
			g_free(line);
			continue;
		}

		/* now read values from each column */
		if (get_value(col[racol], &raval) != SH_SUCCESS) {
			shError("read_data_file: can't read RA value from %s; skipping",
					col[racol]);
			g_free(line);
			continue;
		}
		if (get_value(col[deccol], &decval) != SH_SUCCESS) {
			shError("read_data_file: can't read Dec value from %s; skipping",
					col[deccol]);
			g_free(line);
			continue;
		}

		/*
		 * check to see if the central RA and this star's RA are
		 * wrapped around zero from each other
		 */
		if ((raval < 10) && (central_ra > 350)) {
			delta_ra = (raval + 360) - central_ra;
		} else if ((raval > 350) && (central_ra < 10)) {
			delta_ra = (raval - 360) - central_ra;
		} else {
			delta_ra = raval - central_ra;
		}
		delta_ra *= DEGTORAD;

		/*
		 * let's transform from (delta_RA, delta_Dec) to (xi, eta),
		 */
		dec_rad = decval * DEGTORAD;
		xx = cos(dec_rad) * sin(delta_ra);
		yy = sin(cent_dec_rad) * sin(dec_rad)
				+ cos(cent_dec_rad) * cos(dec_rad) * cos(delta_ra);
		xi = (xx / yy);

		xx = cos(cent_dec_rad) * sin(dec_rad)
				- sin(cent_dec_rad) * cos(dec_rad) * cos(delta_ra);
		eta = (xx / yy);

#ifdef DEBUG
		printf("  xi = %10.5f,   eta = %10.5f\n", xi, eta);
#endif

		/* if desired, convert xi and eta from radians to arcsec */
		if (doASEC > 0) {
			xi *= RADtoASEC;
			eta *= RADtoASEC;
		}
#ifdef DEBUG
		printf("  now  xi = %10.5f,   eta = %10.5f\n", xi, eta);
#endif

		/*
		 * now build up the output line.  Note that we use slightly
		 * different output formats for (xi, eta) if they are expressed
		 * in radians or arcseconds.  
		 */

		line = g_strstrip(line);

		gchar **token = g_strsplit_set(line, " \t", -1);
		i = 0;
		while (token[i] && strcmp(token[i], "\n")) {
			if (i == racol) {
				g_free(token[i]);
				if (doASEC > 0) {
					/* write the offsets in arcsec */
					token[i] = g_strdup_printf("%12.5f", xi);
				} else {
					token[i] = g_strdup_printf("%13.6e", xi);
				}
			} else if (i == deccol) {
				g_free(token[i]);
				if (doASEC > 0) {
					/* write the offsets in arcsec */
					token[i] = g_strdup_printf("%12.5f", eta);
				} else {
					token[i] = g_strdup_printf("%13.6e", eta);
				}
			}
			i++;
		}
		gchar *newline = g_strjoinv(" ", token);
		gchar *output_line = g_strconcat (newline, "\n", NULL);
		GOutputStream *output_stream = (GOutputStream *)g_file_append_to(file_out, G_FILE_CREATE_NONE, NULL, &error);
		if (output_stream) {
			g_output_stream_write_all(output_stream, output_line, strlen(output_line), NULL, NULL, NULL);
		} else {
			if (error != NULL) {
				shError("proc_star_file: can't open file %s for input", file);
				g_clear_error(&error);
			}
		}

		g_object_unref(output_stream);
		g_strfreev(token);
		g_free(newline);
		g_free(output_line);
		g_free(line);
	}

	g_object_unref(data_input);
	g_object_unref(input_stream);
	g_object_unref(file_in);
	return (SH_SUCCESS);
}
