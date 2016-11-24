#ifndef _COMMAND_H_
#define _COMMAND_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

typedef
struct {
	char name[32];
	int nbarg;
	char usage[128];
	int (* process)(int);
} command;

int	process_load(int nb);
#if 0
int	process_convert(int nb);
int	process_trichro(int nb);
int	process_composervb(int nb);
int	process_bmp2fits(int nb);
#endif
int	process_satu(int nb);
int	process_save(int nb);
int	process_savebmp(int nb);

#ifdef HAVE_LIBJPEG
int	process_savejpg(int nb);
#endif
#ifdef HAVE_LIBTIFF
int	process_savetif(int nb);
#endif
int	process_savepnm(int nb);

int	process_imoper(int nb);
int 	process_addmax(int nb);
#if 0
int	process_soper(int nb);
int	process_imoper2(int nb);
int	process_soper2(int nb);
#endif
int	process_fdiv(int nb);
int	process_fmul(int nb);
#if 0
int	process_composit(int nb);
int	process_shift(int nb);
int	process_shift2(int nb);
int	process_rshift2(int nb);
#endif
int	process_entropy(int nb);
int	process_gauss(int nb);
int	process_unsharp(int nb);
#if 0
int	process_unsharp2(int nb);
int	process_fwhm(int nb);
#endif
int	process_crop(int nb);
int	process_seq_crop(int nb);

#if 0
int	process_crop2(int nb);
int	process_medstack(int nb);
int	process_register(int nb);
#endif
int	process_cd(int nb);
int	process_wrecons(int nb);
int	process_wavelet(int nb);
//int	process_animate(int nb);
//
int	process_log(int nb);
int	process_ls(int nb);
int	process_contrast(int nb);
int	process_cdg(int nb);
int 	process_clearstar(int nb);
int	process_mirrorx(int nb);
int	process_mirrory(int nb);
#ifdef HAVE_OPENCV
int	process_resample(int nb);
int	process_rotate(int nb);
#endif
int	process_rotatepi(int nb);
int 	process_psf(int nb);
int	process_seq_psf(int nb);
int	process_bg(int nb);
int	process_bgnoise(int nb);
int	process_histo(int nb);
int	process_thresh(int nb);
int	process_threshlo(int nb);
int	process_threshhi(int nb);
int	process_nozero(int nb);
int	process_ddp(int nb);
int	process_new(int nb);
int	process_visu(int nb);
int	process_fill2(int nb);
int	process_findstar(int nb);
int	process_findhot(int nb);
int	process_cosme(int nb);
int	process_fmedian(int nb);
int	process_fill(int nb);
int	process_offset(int nb);
int	process_scnr(int nb);
int	process_fft(int nb);
int	process_fixbanding(int nb);
int	process_findcosme(int nb);
int	process_split(int nb);
int	process_select(int nb);
int	process_set_mag(int nb);
int	process_set_mag_seq(int nb);
int	process_unset_mag(int nb);
int	process_unset_mag_seq(int nb);
int	process_unselect(int nb);
int	process_stat(int nb);
int	process_stackall(int nb);
#ifdef _OPENMP
int process_set_cpu(int nb);
#endif

#if 0
int	xrgb(int nb, int layer);
int	process_rrgb(int nb);
int	process_grgb(int nb);
int	process_brgb(int nb);
int	process_lrgb(int nb);
#endif
int	process_help(int nb);
int	process_exit(int nb);
int	process_extract(int nb);
int	processcommand(const char *line);

#endif
