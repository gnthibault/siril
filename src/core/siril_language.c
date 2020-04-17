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

/*
 * Language lists that we want to generate only once at program startup:
 * @l10n_lang_list: all available localizations self-localized;
 * @all_lang_list: all known languages, in the user-selected language.
 */

#include <locale.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/siril_app_dirs.h"
#include "algos/sorting.h"
#include "gui/callbacks.h"

#include "siril_language.h"

static GHashTable *l10n_lang_list = NULL;
static GHashTable *string_lang_list = NULL;

parsed_code locale_str_test[] = {
	{"de", "Deutsch"},
	{"el", "Greek"},
	{"en", "English"},
	{"es_ES", "Espanol"},
	{"fr", "Français"},
	{"it_IT", "Italiano"},
	{"nl_BE", "Nederlands"},
	{"pl_PL", "Polish"},
	{"ru", "Russian"},
	{"zh_CN", "Chinese (Simplified)"},
	{NULL, NULL}
};

static GHashTable *parse_locale_codes(GHashTable *table) {
	GHashTableIter lang_iter;
	gpointer key;

	GHashTable *string_lang = g_hash_table_new_full(g_str_hash, g_str_equal,
			(GDestroyNotify) g_free, (GDestroyNotify) g_free);

	g_hash_table_iter_init(&lang_iter, table);
	while (g_hash_table_iter_next(&lang_iter, &key, NULL)) {
		gchar *code = (gchar*) key;
		gchar *str_name = NULL;
		int i = 0;

		while (locale_str_test[i].locale) {
			if (!g_strcmp0(code, locale_str_test[i].locale)) {
				str_name = g_strdup(locale_str_test[i].language_name);
				break;
			}
			i++;
		}
		g_hash_table_insert(string_lang,
				g_strdup_printf("%s [%s]", str_name ? str_name : "???", code),
				NULL);
	}
	return string_lang;
}

/* extract locale from a string following this pattern:
 * xxxxxxxxx [locale]
 */
static gchar *extract_locale_from_string(gchar *str) {
	gchar *locale = g_strstr_len(str, -1, "[");
	return g_strndup(locale + 1, strlen(locale) - 2);
}

void siril_language_parser_init() {
	GDir *locales_dir;

	if (l10n_lang_list != NULL) {
		g_warning("siril_language_parser_init() must be run only once.");
		return;
	}

	l10n_lang_list = g_hash_table_new_full(g_str_hash, g_str_equal,
			(GDestroyNotify) g_free, (GDestroyNotify) g_free);

	/* Check all locales we have translations for. */
	locales_dir = g_dir_open(siril_get_locale_dir(), 0, NULL);
	if (locales_dir) {
		const gchar *locale;

		while ((locale = g_dir_read_name(locales_dir)) != NULL) {
			gchar *filename = g_build_filename(siril_get_locale_dir(), locale,
					"LC_MESSAGES",
					PACKAGE ".mo",
					NULL);
			if (g_file_test(filename, G_FILE_TEST_EXISTS)) {
				/* Save the full language code. */
				g_hash_table_insert(l10n_lang_list, g_strdup(locale), NULL);
			}

			g_free(filename);
		}

		g_dir_close(locales_dir);
	}
	string_lang_list = parse_locale_codes(l10n_lang_list);
}

void siril_language_fill_combo() {
	GtkComboBoxText *lang_combo = GTK_COMBO_BOX_TEXT(lookup_widget("combo_language"));
	GList *list = g_hash_table_get_keys(string_lang_list);
	gboolean lang_changed = FALSE;
	int i = 1;

	list = g_list_sort(list, (GCompareFunc) strcompare);

	for (GList *l = list; l; l = l->next) {
		gtk_combo_box_text_append_text(lang_combo, l->data);
		gchar *locale = extract_locale_from_string(l->data);
		if (!g_strcmp0(com.combo_lang, locale)) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(lang_combo), i);
			lang_changed = TRUE;
		}
		g_free(locale);
		i++;
	}
	if (!lang_changed) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(lang_combo), 0);
	}
}

void language_init(const gchar *language) {
	if (!language)
		return;

	g_setenv("LANGUAGE", language, TRUE);
	setlocale(LC_ALL, "");
}

void update_language() {
	GtkComboBoxText *lang_combo = GTK_COMBO_BOX_TEXT(lookup_widget("combo_language"));

	g_free(com.combo_lang);
	if (gtk_combo_box_get_active(GTK_COMBO_BOX(lang_combo)) == 0) {
		com.combo_lang = g_strdup("");
	} else {
		gchar *str = gtk_combo_box_text_get_active_text(lang_combo);
		com.combo_lang = extract_locale_from_string(str);
	}
	writeinitfile();
}
