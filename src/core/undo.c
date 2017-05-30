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

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include "core/siril.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "core/undo.h"
#include "gui/histogram.h"
#include "core/proto.h"

/* *filename must be freed */
static int undo_build_swapfile(fits *fit, char **filename) {
	char *nameBuff;
	char name[] = "/siril_swp-XXXXXX";
	char *tmpdir;
	int len, fd, size;

	tmpdir = com.swap_dir;
	len = strlen(tmpdir) + strlen(name) + 1;

	nameBuff = calloc(1, len * sizeof(char));

	snprintf(nameBuff, len, "%s%s", tmpdir, name);
	fd = mkstemp(nameBuff);
	if (fd < 1) {
		siril_log_message(_("File I/O Error: Unable to create swap file in %s: [%s]\n"),
				tmpdir, strerror(errno));
		free(nameBuff);
		return 1;
	}

	size = fit->rx * fit->ry * fit->naxes[2];

	errno = 0;
	// Write some data to the temporary file
	if (-1 == write(fd, fit->data, size * sizeof(WORD))) {
		siril_log_message(_("File I/O Error: Unable to write swap file in %s: [%s]\n"),
				tmpdir, strerror(errno));
		free(nameBuff);
		close(fd);
		return 1;
	}
	*filename = nameBuff;
	close(fd);

	return 0;
}

static int undo_remove_item(historic *histo, int index) {
	if (histo[index].filename) {
		unlink(histo[index].filename);
		free(histo[index].filename);
		histo[index].filename = NULL;
	}
	memset(histo[index].history, 0, FLEN_VALUE);
	return 0;
}

static void undo_add_item(fits *fit, char *filename, char *histo) {

	if (!com.history) {
		com.hist_size = HISTORY_SIZE;
		com.history = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}
	/* when undo, we remove all further items being after */
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}
	com.history[com.hist_current].filename = filename;
	com.history[com.hist_current].rx = fit->rx;
	com.history[com.hist_current].ry = fit->ry;
	snprintf(com.history[com.hist_current].history, FLEN_VALUE, "%s", histo);

	if (com.hist_current == com.hist_size - 1) {
		/* we must shift all elements except 0 that must always match with the original file
		 * 0  1  2  3  4  5  6  7  8  9 10 become
		 * 0  2  3  4  5  6  7  8  9 10 11 and
		 * 0  3  4  5  6  7  8  9 10 11 12 and so on
		 */
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
				(com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}
	com.hist_current++;
	com.hist_display = com.hist_current;
}

static int undo_get_data(fits *fit, historic hist) {
	int fd;
	unsigned int size;
	WORD *buf;

	if ((fd = open(hist.filename, O_RDONLY)) == -1) {
		printf("Error opening swap file : %s\n", hist.filename);
		return 1;
	}

	errno = 0;
	fit->rx = hist.rx;
	fit->ry = hist.ry;
	size = fit->rx * fit->ry * fit->naxes[2];
	buf = calloc(1, size * sizeof(WORD));
	// read the data from temporary file
	if ((read(fd, buf, size * sizeof(WORD)) < size * sizeof(WORD))) {
		printf("Read failed with error [%s]\n", strerror(errno));
		free(buf);
		close(fd);
		return 1;
	}
	WORD *newdata = (WORD*) realloc(fit->data, size * sizeof(WORD));
	if (!newdata) {
		free(newdata);
		free(buf);
		close(fd);
		return 1;
	}
	fit->data = newdata;
	memcpy(fit->data, buf, size * sizeof(WORD));
	fit->pdata[RLAYER] = fit->data;
	if (fit->naxes[2] > 1) {
		fit->pdata[GLAYER] = fit->data + fit->rx * fit->ry;
		fit->pdata[BLAYER] = fit->data + fit->rx * fit->ry * 2;
	}
	free(buf);
	close(fd);
	return 0;
}

gboolean is_undo_available() {
    return (com.history && com.hist_display > 0);
}

gboolean is_redo_available() {
    return (com.history && (com.hist_display < com.hist_current - 1));
}

int undo_save_state(char *message, ...) {
	char *filename;
	char histo[FLEN_VALUE];
	va_list args;
	va_start(args, message);

	if (single_image_is_loaded()) {
		if (message == NULL)
			memset(histo, 0, FLEN_VALUE);
		else
			vsnprintf(histo, FLEN_VALUE, message, args);

		if (undo_build_swapfile(&gfit, &filename)) {
			va_end(args);
			return 1;
		}

		undo_add_item(&gfit, filename, histo);

		/* update menus */
		update_MenuItem();
	}
	va_end(args);
	return 0;
}

int undo_display_data(int dir) {
	if (!com.history) {
		return 1;
	}
	switch (dir) {
	case UNDO:
		if (is_undo_available()) {
			if (com.hist_current == com.hist_display) {
				undo_save_state(NULL);
				com.hist_display--;
			}
			com.hist_display--;
			undo_get_data(&gfit, com.history[com.hist_display]);
			update_gfit_histogram_if_needed();
			redraw(com.cvport, REMAP_ALL);
			redraw_previews();
		}
		break;
	case REDO:
		if (is_redo_available()) {
			com.hist_display++;
			undo_get_data(&gfit, com.history[com.hist_display]);
			update_gfit_histogram_if_needed();
			redraw(com.cvport, REMAP_ALL);
			redraw_previews();
		}
		break;
	default:
		printf("ERROR\n");
		return -1;
	}
	return 0;
}

int undo_flush() {
	int i;

	if (!com.history) {
		return 1;
	}
	for (i = 0; i < com.hist_current; i++) {
		undo_remove_item(com.history, i);
	}
	free(com.history);
	com.history = NULL;
	com.hist_current = 0;
	com.hist_display = 0;
	return 0;
}
