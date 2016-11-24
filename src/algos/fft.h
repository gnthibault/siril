#ifndef FFT_H
#define FFT_H

/* fft data from GUI */
struct fft_data {
	fits *fit;
	char *type;
	const char *modulus, *phase;
	int type_order;
};

void fft_to_spectra(fits *, fftw_complex *, double *, double *);
void fft_to_fr(fits *, fftw_complex *, double *, double *);
void change_symmetry(fits *, unsigned int, unsigned int, unsigned int *,
		unsigned int *);
double normalisation_spectra(fits *, double *, double *, WORD *, WORD *, int);
void save_dft_information_in_gfit(fits *);
void FFTD(fits *, fits *, fits *, int, int);
void FFTI(fits*, fits*, fits*, int, int);
gpointer fourier_transform(gpointer p);
#endif
