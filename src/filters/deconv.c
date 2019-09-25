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


#include "core/siril.h"
#include "core/proto.h"

#include "core/undo.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "opencv/opencv.h"

#include "deconv.h"

gpointer LRdeconv(gpointer p) {
	struct RL_data *args = (struct RL_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Lucy-Richardson deconvolution: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	cvLucyRichardson(args->fit, args->sigma, args->iter);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	siril_add_idle(end_generic, args);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

/************************ GUI for deconvolution ***********************/
void on_menuitem_deconvolution_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("deconvolution_dialog");
}

void on_deconvolution_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("deconvolution_dialog");
}

void on_deconvolution_Apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *iter;
	GtkRange *sigma;

	iter = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"spin_iterations"));
	sigma = GTK_RANGE(gtk_builder_get_object(builder, "scale_deconvolution"));

	struct RL_data *args = malloc(sizeof(struct RL_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->sigma = gtk_range_get_value(sigma);
	args->iter = gtk_spin_button_get_value(iter);

	undo_save_state(&gfit, "Processing: Deconvolution (iter=%d, sig=%.3f)", args->iter,
			args->sigma);
	start_in_new_thread(LRdeconv, args);
}
