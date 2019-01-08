#include "../core/siril.h"
#include "../core/proto.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
	fits fits1 = {0}, fits2 = {0};

	if (argc < 3) {
		fprintf(stderr, "Usage: %s image1.fit image2.fit\n", *argv);
		exit(2);
	}

	if (readfits(argv[1], &fits1, NULL) || readfits(argv[2], &fits2, NULL)) {
		exit(2);
	}

	if (fits1.header && fits2.header) {
		if (strcmp(fits1.header, fits2.header))
			fprintf(stdout, "headers differ\n");
	}
	if (fits1.naxis != fits2.naxis) {
		fprintf(stdout, "number of axis differ\n");
		exit(1);
	}
	if (fits1.naxes[0] != fits2.naxes[0] || fits1.naxes[1] != fits2.naxes[1] ||
			fits1.naxes[2] != fits2.naxes[2]) {
		fprintf(stdout, "image axis differ\n");
		exit(1);
	}

	if (memcmp(fits1.data, fits2.data, fits1.naxes[0] * fits1.naxes[1] * fits1.naxes[2])) {
		fprintf(stdout, "image data differ\n");
		exit(1);
	}

	fprintf(stdout, "images are identical\n");
	return 0;
}
