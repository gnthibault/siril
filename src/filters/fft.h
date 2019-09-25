#ifndef FFT_H
#define FFT_H

/* fft data from GUI */
struct fft_data {
	fits *fit;
	char *type;
	const char *modulus, *phase;
	int type_order;
	int retval;
};

gpointer fourier_transform(gpointer p);
#endif
