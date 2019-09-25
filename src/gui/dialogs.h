#ifndef SRC_GUI_DIALOGS_H_
#define SRC_GUI_DIALOGS_H_

#include "core/siril.h"
#include "core/proto.h"

typedef enum {
	NO_DIALOG = -1,
	INFORMATION_DIALOG,
	IMAGE_PROCESSING_DIALOG,
	OTHER_DIALOG
} DialogType;

struct _SirilDialogEntry {
	gchar *identifier;
	DialogType type;

	gboolean has_preview;

	void (*apply_function)(void);
};

void siril_open_dialog(gchar *id);
void siril_close_dialog(gchar *id);
void siril_close_preview_dialogs();

#endif /* SRC_GUI_DIALOGS_H_ */
