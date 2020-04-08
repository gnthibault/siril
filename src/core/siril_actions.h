/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_CORE_SIRIL_ACTIONS_H_
#define SRC_CORE_SIRIL_ACTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <glib.h>

void open_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void cwd_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void save_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void save_as_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void undo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void redo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void quit_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void about_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void preferences_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void close_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void scripts_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) ;
#ifdef HAVE_LIBCURL
void updates_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
#endif
void full_screen_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void keyboard_shortcuts_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_conversion_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_sequence_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_prepro_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_registration_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_plot_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_stacking_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) ;
void tab_logs_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void toolbar_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) ;

#endif /* SRC_CORE_SIRIL_ACTIONS_H_ */
