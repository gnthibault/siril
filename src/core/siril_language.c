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
#include "gui/callbacks.h"

#include "siril_language.h"

static GHashTable *l10n_lang_list = NULL;
static GHashTable *full_lang_list = NULL;

parsed_code locale_str[] = {
	{"ar_DZ", "العربية"},
	{"de", "Deutsch"},
	{"el", "Ελληνικά"},
	{"en", "English"},
	{"es_ES", "Espanol"},
	{"fr", "Français"},
	{"it_IT", "Italiano"},
	{"ja_JP", "日本語"},
	{"nl_BE", "Nederlands"},
	{"pl_PL", "Polish"},
	{"pt_PT", "Português"},
	{"ru", "русский"},
	{"zh_CN", "汉语"},
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

		while (locale_str[i].locale) {
			if (!g_strcmp0(code, locale_str[i].locale)) {
				str_name = g_strdup(locale_str[i].language_name);
				break;
			}
			i++;
		}
		g_hash_table_insert(string_lang,
				g_strdup_printf("%s [%s]", str_name ? str_name : "???", code),
				NULL);
		g_free(str_name);
	}
	return string_lang;
}

/* extract locale from a string following this pattern:
 * xxxxxxxxx [locale]
 */
static gchar *extract_locale_from_string(const gchar *str) {
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

	/* By default Siril is written in english */
	g_hash_table_insert(l10n_lang_list, g_strdup("en"), NULL);
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
	full_lang_list = parse_locale_codes(l10n_lang_list);
}

static gint locale_compare(gconstpointer *a, gconstpointer *b) {
	const gchar *s1 = (const gchar *) a;
	const gchar *s2 = (const gchar *) b;

	gchar *k1 = extract_locale_from_string(s1);
	gchar *k2 = extract_locale_from_string(s2);

	gint result = g_strcmp0(k1, k2);
	g_free(k1);
	g_free(k2);

	return result;
}

void siril_language_fill_combo(const gchar *language) {
	GtkComboBoxText *lang_combo = GTK_COMBO_BOX_TEXT(lookup_widget("combo_language"));
	GList *list = g_hash_table_get_keys(full_lang_list);
	gboolean lang_changed = FALSE;
	int i = 1;

	gtk_combo_box_text_remove_all(lang_combo);
	gtk_combo_box_text_append(lang_combo, 0, _("System Language"));

	list = g_list_sort(list, (GCompareFunc) locale_compare);

	for (GList *l = list; l; l = l->next) {
		gtk_combo_box_text_append_text(lang_combo, l->data);
		gchar *locale = extract_locale_from_string(l->data);
		if (!g_strcmp0(language, locale)) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(lang_combo), i);
			lang_changed = TRUE;
		}
		g_free(locale);
		i++;
	}
	if (!lang_changed) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(lang_combo), 0);
	}

	g_list_free(list);
}

void language_init(const gchar *language) {
	if ((!language) || (language[0] == '\0'))
		return;

	/* This is default language */
	if (!g_ascii_strcasecmp(language, "en")) {
		g_setenv("LANGUAGE", "C", TRUE);
	} else {
		g_setenv("LANGUAGE", language, TRUE);
	}
	setlocale(LC_ALL, "");
}

gchar *get_interface_language() {
	GtkComboBoxText *lang_combo = GTK_COMBO_BOX_TEXT(lookup_widget("combo_language"));

	if (gtk_combo_box_get_active(GTK_COMBO_BOX(lang_combo)) == 0) {
		return g_strdup("");
	} else {
		gchar *str = gtk_combo_box_text_get_active_text(lang_combo);
		return extract_locale_from_string(str);
	}
}
