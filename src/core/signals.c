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

#define ANSI_COLOR_RED     "\e[1m\x1b[31m"
#define ANSI_COLOR_RESET   "\x1b[0m\e[0m"

#include <signal.h>
#ifdef _WIN32
#include <windows.h>
#include <dbghelp.h>
#else
#include <execinfo.h>
#endif

#include "siril.h"
#include "proto.h"
#include "signals.h"

#define STACK_DEPTH 256

static void signal_handled(int s) {
	gchar *visit = NULL;

	g_printf("Error, signal %d:\n", s);
	visit = g_strdup_printf(_("Please report this bug to: %s"), PACKAGE_BUGREPORT);
	switch (s) {
	case SIGSEGV:
	case SIGFPE:
	case SIGABRT:
	case SIGILL:
		g_printf(ANSI_COLOR_RED"%s\n"ANSI_COLOR_RESET, visit);
	}

#if (!defined _WIN32 && defined HAVE_EXECINFO_H)
	int i;
	void *stack[STACK_DEPTH];

	size_t size = backtrace(stack, sizeof(stack) / sizeof(void*));

	char **message = backtrace_symbols(stack, size);
	if (message != NULL && message[0] != NULL) {
		for (i = 0; i < size && message != NULL; ++i) {
			g_printf("[#%i] in %s\n", i, message[i]);
		}
		free(message);
	}
#else
	unsigned int i;
	void *stack[STACK_DEPTH];
	unsigned short size;
	SYMBOL_INFO *symbol;
	HANDLE process;

	process = GetCurrentProcess();

	SymInitialize(process, NULL, TRUE);

	size = CaptureStackBackTrace(0, sizeof(stack) / sizeof(void*), stack, NULL);
	symbol = (SYMBOL_INFO*) calloc(sizeof(SYMBOL_INFO) + 256 * sizeof(char), 1);
	symbol->MaxNameLen = 255;
	symbol->SizeOfStruct = sizeof(SYMBOL_INFO);

	for (i = 0; i < size; i++) {
		SymFromAddr(process, (DWORD64)(stack[i]), 0, symbol);

		g_printf("[#%i]: in %s\n", i, symbol->Name);
	}

	free(symbol);
#endif
	g_free(visit);
	gtk_main_quit();
}

void signals_init() {
#ifndef _WIN32
	signal(SIGHUP, signal_handled);
	signal(SIGQUIT, signal_handled);
	signal(SIGBUS, signal_handled);
	signal(SIGINT, signal_handled);
	signal(SIGTRAP, signal_handled);
#endif
	signal(SIGABRT, signal_handled);
	signal(SIGFPE, signal_handled);
	signal(SIGSEGV, signal_handled);
	signal(SIGTERM, signal_handled);
	signal(SIGILL, signal_handled);
}
