#ifndef SRC_GUI_SCRIPT_MENU_H_
#define SRC_GUI_SCRIPT_MENU_H_

GSList *get_list_from_preferences_dialog();
GSList *set_list_to_preferences_dialog(GSList *list);
int initialize_script_menu();
int refresh_scripts(gboolean update_list, gchar **error);
void siril_get_on_script_pages();

#endif /* SRC_GUI_SCRIPT_MENU_H_ */
