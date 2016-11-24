#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

void get_stars_list_store();
void add_star_to_list(fitted_PSF *);
void fill_stars_list(fitted_PSF **);
void refresh_stars_list(fitted_PSF **);
void clear_stars_list();
void display_PSF(fitted_PSF **);
void remove_selected_line();
void move_selected_line();
void remove_all_lines();
void display_status();

#endif
