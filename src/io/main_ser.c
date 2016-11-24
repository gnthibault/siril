#include "io/ser.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
	struct ser_struct ser_file;
	ser_init_struct(&ser_file);
	if (argc == 1) {
		fprintf(stdout, "Usage: %s ser_file_name\n", *argv);
		exit(1);
	}
	ser_open_file(argv[1], &ser_file);
	display_ser_info(&ser_file);
	return 0;
}
