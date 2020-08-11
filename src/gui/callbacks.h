#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <sys/time.h>
#include <stdint.h>
#include "core/siril.h"	// for sliders_mode

GtkWidget* lookup_widget (const gchar *widget_name);

void initialize_all_GUI(gchar *files);
void load_prefered_theme();
void set_cutoff_sliders_max_values();		// was set_upper_minmax
void set_cutoff_sliders_values();		// was set_ranges
void set_sliders_value_to_gfit();
void initialize_display_mode();
void set_display_mode();
void adjust_exclude(int n, gboolean changed);
int adjust_sellabel();
void set_GUI_CWD();
void set_GUI_MEM(unsigned long long size);
void set_GUI_DiskSpace(int64_t mem);
void set_GUI_misc();
void set_icon_entry(GtkEntry *entry, gchar *string);
void update_MenuItem();
void sliders_mode_set_state(sliders_mode);
int copy_rendering_settings_when_chained(gboolean from_GUI);

void clear_sampling_setting_box();
void set_GUI_CAMERA();

int match_drawing_area_widget(GtkWidget *drawing_area, gboolean allow_rgb);
const char *vport_number_to_name(int);
const char *untranslated_vport_number_to_name(int);
void update_display_selection();
void update_display_fwhm();
void display_filename();
void set_precision_switch();
void set_layers_for_assign();
void set_layers_for_registration();
void show_dialog(const char *text, const char *title, const char *icon);
void show_txt_and_data_dialog(const char *text, const char *data, const char *title, const char *icon);
void show_data_dialog(char *text, char *title);
GtkWindow *siril_get_active_window();
void initialize_FITS_name_entries();

void adjust_vport_size_to_image();
void scrollbars_hadjustment_changed_handler(GtkAdjustment *adjustment, gpointer user_data);
void scrollbars_vadjustment_changed_handler(GtkAdjustment *adjustment, gpointer user_data);
void set_output_filename_to_sequence_name();
void close_tab();
void activate_tab(int vport);
void control_window_switch_to_tab(main_tabs tab);

void update_prepro_interface();

void update_statusbar_convert();

void update_spinCPU(int max);

void save_main_window_state();
void load_main_window_state();
void siril_quit();

/* for image_display */
void set_viewer_mode_widgets_sensitive(gboolean sensitive);

/*****************************************************************************
*      P U B L I C      C A L L B A C K      F U N C T I O N S               *
 ****************************************************************************/
void on_radiobutton_minmax_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_radiobutton_hilo_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_radiobutton_user_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_checkchain_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_max_entry_changed(GtkEditable *editable, gpointer user_data);
void on_min_entry_changed(GtkEditable *editable, gpointer user_data);

void on_seqproc_entry_changed (GtkComboBox *widget,	gpointer user_data);
void on_excludebutton_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data);

void on_spin_w_changed(GtkSpinButton *spinbutton, gpointer user_data);

void on_checkbutton_cam_toggled(GtkButton *button, gpointer user_data);
void on_checkbutton_auto_toggled(GtkButton *button, gpointer user_data);

gboolean on_drawingarea_key_press_event(GtkWidget *widget, GdkEventKey *event, gpointer user_data);

#endif
