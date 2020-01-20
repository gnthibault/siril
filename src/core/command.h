#ifndef SRC_CORE_COMMAND_H_
#define SRC_CORE_COMMAND_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

int	process_addmax(int nb);
int	process_asinh(int nb);

int	process_bg(int nb);
int	process_bgnoise(int nb);

int	process_cd(int nb);
int	process_cdg(int nb);
int	process_crop(int nb);
int	process_clahe(int nb);
int	process_clear(int nb);
int	process_clearstar(int nb);
int	process_close(int nb);
int	process_convertraw(int nb);
int	process_cosme(int nb);

int	process_ddp(int nb);

int	process_entropy(int nb);
int	process_exit(int nb);
int	process_extract(int nb);

int	process_fdiv(int nb);
int	process_fft(int nb);
int	process_fill(int nb);
int	process_fill2(int nb);
int	process_findcosme(int nb);
int	process_findhot(int nb);
int	process_findstar(int nb);
int	process_fixbanding(int nb);
int	process_fmedian(int nb);
int	process_fmul(int nb);

int	process_gauss(int nb);
int	process_grey_flat(int nb);

int	process_help(int nb);
int	process_histo(int nb);

int	process_imoper(int nb);

int	process_load(int nb);
int	process_log(int nb);
int	process_ls(int nb);

int	process_mirrorx(int nb);
int	process_mirrory(int nb);

int	process_new(int nb);
int	process_nozero(int nb);

int	process_offset(int nb);

int	process_preprocess(int nb);
int	process_psf(int nb);

int	process_register(int nb);
int	process_resample(int nb);
int	process_rgradient(int nb);
int	process_rl(int nb);
int	process_rotate(int nb);
int	process_rotatepi(int nb);

int	process_satu(int nb);
int	process_save(int nb);
int	process_savebmp(int nb);
#ifdef HAVE_LIBJPEG
int	process_savejpg(int nb);
#endif
#ifdef HAVE_LIBPNG
int	process_savepng(int nb);
#endif
int	process_savepnm(int nb);
#ifdef HAVE_LIBTIFF
int	process_savetif(int nb);
#endif
int	process_scnr(int nb);
int	process_select(int nb);
#ifdef _OPENMP
int	process_set_cpu(int nb);
#endif
int	process_set_ext(int nb);
int	process_set_findstar(int nb);
int	process_set_mag(int nb);
int	process_set_mag_seq(int nb);
int	process_set_mem(int nb);
int process_set_ref(int nb);
int	process_seq_crop(int nb);
int	process_seq_psf(int nb);
int	process_seq_split_cfa(int nb);
int	process_split(int nb);
int	process_split_cfa(int nb);
int	process_stat(int nb);
int	process_stackall(int nb);
int	process_stackone(int nb);

int	process_thresh(int nb);
int	process_threshlo(int nb);
int	process_threshhi(int nb);

int	process_unset_mag(int nb);
int	process_unset_mag_seq(int nb);
int	process_unselect(int nb);
int	process_unsharp(int nb);

int	process_visu(int nb);

int	process_wavelet(int nb);
int	process_wrecons(int nb);

#endif
