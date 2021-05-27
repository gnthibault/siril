#ifndef SRC_CORE_COMMAND_LIST_H_
#define SRC_CORE_COMMAND_LIST_H_


#include "siril.h"
#include "command.h"
#include "command_def.h"

#define MAX_COMMAND_WORDS 50		// max number of words to split in command line input

extern char *word[MAX_COMMAND_WORDS];	// NULL terminated

typedef
struct {
	char *name;
	int nbarg;
	char *usage;
	int (* process)(int);
	char *definition;
	gboolean scriptable;
} command;

static command commands[] = {
	/* name,	nbarg,	usage,		function pointer, definition, scriptable */
	{"addmax", 1,	"addmax filename", process_addmax, STR_ADDMAX, FALSE},
	{"asinh", 1,	"asinh stretch", process_asinh, STR_ASINH, TRUE},

	{"bg", 0, "bg", process_bg, STR_BG, TRUE},
	{"bgnoise", 0, "bgnoise", process_bgnoise, STR_BGNOISE, TRUE},

	{"cd", 1, "cd directory", process_cd, STR_CD, TRUE},
	{"cdg", 0, "cdg", process_cdg, STR_CDG, TRUE},
	{"clahe", 2, "clahe cliplimit tileSize", process_clahe, STR_CLAHE, TRUE},
	{"clear", 0, "clear", process_clear, STR_CLEAR, FALSE},
	{"clearstar", 0, "clearstar", process_clearstar, STR_CLEARSTAR, FALSE},
	{"close", 0, "close", process_close, STR_CLOSE, TRUE},
	{"convert", 1, "convert basename [-debayer] [-fitseq] [-start=index] [-out=]", process_convert, STR_CONVERT, TRUE},
	{"convertraw", 1, "convertraw basename [-debayer] [-fitseq] [-start=index] [-out=]", process_convertraw, STR_CONVERTRAW, TRUE},
	{"cosme", 1, "cosme [filename].lst", process_cosme, STR_COSME, TRUE},
	{"cosme_cfa", 1, "cosme_cfa [filename].lst", process_cosme, STR_COSME_CFA, TRUE},
	{"crop", 0, "crop [x y width height]", process_crop, STR_CROP, TRUE},

	{"ddp", 3, "ddp level coef sigma", process_ddp, STR_DDP, FALSE},

	{"entropy", 0, "entropy", process_entropy, STR_ENTROPY, TRUE},
	{"exit", 0, "exit", process_exit, STR_EXIT, TRUE},
	{"extract", 1, "extract NbPlans", process_extract, STR_EXTRACT, TRUE},
	{"extract_Ha", 0, "extract_Ha", process_extractHa, STR_EXTRACTHA, TRUE},
	{"extract_HaOIII", 0, "extract_HaOIII", process_extractHaOIII, STR_EXTRACTHAOIII, TRUE},

	{"fdiv", 2, "fdiv filename scalar", process_fdiv, STR_FDIV, TRUE},
	{"fftd", 2, "fftd modulus phase", process_fft, STR_FFTD, TRUE},
	{"ffti", 2, "ffti modulus phase", process_fft, STR_FFTI, TRUE},
	{"fill", 1, "fill value [x y width height]", process_fill, STR_FILL, TRUE},
	{"fill2", 1, "fill2 value [x y width height]", process_fill2, STR_FILL2, TRUE},
	{"find_cosme", 2, "find_cosme cold_sigma hot_sigma", process_findcosme, STR_FIND_COSME, TRUE},
	{"find_cosme_cfa", 2, "find_cosme_cfa cold_sigma hot_sigma", process_findcosme, STR_FIND_COSME_CFA, TRUE},
	{"find_hot", 3, "find_hot filename cold_sigma hot_sigma", process_findhot, STR_FIND_HOT, TRUE},
	{"findstar", 0, "findstar", process_findstar, STR_FINDSTAR, FALSE},
	{"fix_xtrans", 0, "fix_xtrans", process_fix_xtrans, STR_FIXXTRANS, TRUE},
	{"fixbanding", 2, "fixbanding amount sigma", process_fixbanding, STR_FIXBANDING, TRUE},
	{"fmedian", 2, "fmedian ksize modulation", process_fmedian, STR_FMEDIAN, TRUE},
	{"fmul", 1, "fmul scalar", process_fmul, STR_FMUL, TRUE},

	{"gauss", 1, "gauss sigma", process_gauss, STR_GAUSS, TRUE},
	{"grey_flat", 0, "grey_flat", process_grey_flat, STR_GREY_FLAT, TRUE},

	{"help", 0, "help", process_help, STR_HELP, FALSE},
	{"histo", 1, "histo channel (channel=0, 1, 2 with 0: red, 1: green, 2: blue)", process_histo, STR_HISTO, TRUE},

	/* commands oper filename and curent image */
	{"iadd", 1, "iadd filename", process_imoper, STR_IADD, TRUE},
	{"idiv", 1, "idiv filename", process_imoper, STR_IDIV, TRUE},
	{"imul", 1, "imul filename", process_imoper, STR_IMUL, TRUE},
	{"isub", 1, "isub filename", process_imoper, STR_ISUB, TRUE},

	{"linear_match", 2, "linear_match reference low high", process_linear_match, STR_LMATCH, TRUE}, /* logarifies current image */
	{"link", 1, "link basename [-start=index] [-out=]", process_link, STR_LINK, TRUE},
	{"load", 1, "load filename.[ext]", process_load, STR_LOAD, TRUE},
	// specific loads are not required, but could be used to force the
	// extension to a higher priority in case two files with same basename
	// exist (stat_file() manages that priority order for now).
	{"log", 0, "log", process_log, STR_LOG, TRUE}, /* logarifies current image */
#ifndef _WIN32
	{"ls", 0, "ls", process_ls, STR_LS, FALSE},
#endif

	{"merge", 3, "merge sequence1 sequence2 [sequence3 ...] output_sequence", process_merge, STR_MERGE, TRUE},
	{"mirrorx", 0, "mirrorx", process_mirrorx, STR_MIRRORX, TRUE},
	{"mirrory", 0, "mirrory", process_mirrory, STR_MIRRORY, TRUE},
	{"mtf", 3, "mtf low mid high", process_mtf, STR_MTF, TRUE},

	{"neg", 0, "neg", process_neg, STR_NEG, TRUE},
	{"new", 3, "new width height nb_channel", process_new, STR_NEW, FALSE},
	{"nozero", 1, "nozero level (replaces null values by level)", process_nozero, STR_NOZERO, TRUE}, /* replaces null values by level */

	{"offset", 1, "offset value", process_offset, STR_OFFSET, TRUE},

	{"preprocess", 1, "preprocess sequencename [-bias=filename] [-dark=filename] [-flat=filename] [-cfa] [-debayer] [-flip] [-equalize_cfa] [-opt] [-prefix=]", process_preprocess, STR_PREPROCESS, TRUE},
	{"psf", 0, "psf", process_psf, STR_PSF, FALSE},

	{"register", 1, "register sequence [-norot] [-drizzle] [-prefix=] [-minpairs=]", process_register, STR_REGISTER, TRUE},
	{"reloadscripts", 0, "reloadscripts", process_reloadscripts, STR_RELOADSCRIPTS, FALSE},
	{"requires", 1, "requires", process_requires, STR_REQUIRES, TRUE},
	{"resample", 1, "resample factor", process_resample, STR_RESAMPLE, TRUE},
	{"rgradient", 4, "rgradient xc yc dR dalpha", process_rgradient, STR_RGRADIENT, TRUE},
	{"rl", 3, "rl sigma corner_radius_boost iterations", process_rl, STR_RL, TRUE},
	{"rmgreen", 1, "rmgreen type", process_scnr, STR_RMGREEN, TRUE},
	{"rotate", 1, "rotate degree [-nocrop]", process_rotate, STR_ROTATE, TRUE},
	{"rotatePi", 0, "rotatePi", process_rotatepi, STR_ROTATEPI, TRUE},

	{"satu", 1, "satu coeff", process_satu, STR_SATU, TRUE},
	{"save", 1, "save filename", process_save, STR_SAVE, TRUE},
	{"savebmp", 1, "savebmp filename", process_savebmp, STR_SAVEBMP, TRUE},
#ifdef HAVE_LIBJPEG
	{"savejpg", 1, "savejpg filename [quality]", process_savejpg, STR_SAVEJPG, TRUE},
#endif
#ifdef HAVE_LIBPNG
	{"savepng", 1, "savepng filename", process_savepng, STR_SAVEPNG, TRUE},
#endif
	{"savepnm", 1, "savepnm filename", process_savepnm, STR_SAVEPNM, TRUE},
#ifdef HAVE_LIBTIFF
	{"savetif", 1, "savetif filename", process_savetif, STR_SAVETIF, TRUE},
	{"savetif32", 1, "savetif32 filename", process_savetif, STR_SAVETIF32, TRUE},
	{"savetif8", 1, "savetif8 filename", process_savetif, STR_SAVETIF8, TRUE},
#endif
	{"select", 2, "select from to", process_select, STR_SELECT, FALSE},
	{"seqextract_Ha", 1, "seqextractHa sequencename [-prefix=]", process_seq_extractHa, STR_SEQEXTRACTHA, TRUE},
	{"seqextract_HaOIII", 1, "seqextractHaOIII sequencename", process_seq_extractHaOIII, STR_SEQEXTRACTHAOIII, TRUE},
	{"seqcrop", 0, "seqcrop [x y width height] [-prefix=]", process_seq_crop, STR_SEQCROP, FALSE},
	{"seqfind_cosme", 3, "seqfind_cosme sequencename cold_sigma hot_sigma [-prefix=]", process_findcosme, STR_SEQFIND_COSME, TRUE},
	{"seqfind_cosme_cfa", 3, "seqfind_cosme_cfa sequencename cold_sigma hot_sigma [-prefix=]", process_findcosme, STR_SEQFIND_COSME_CFA, TRUE},
	{"seqmtf", 4, "seqmtf sequencename low mid high [-prefix=]", process_seq_mtf, STR_SEQMTF, TRUE},
	{"seqpsf", 0, "seqpsf", process_seq_psf, STR_SEQPSF, FALSE},
	{"seqsplit_cfa", 1, "seqsplit_cfa sequencename [-prefix=]", process_seq_split_cfa, STR_SEQSPLIT_CFA, TRUE},
	{"seqstat", 2, "seqstat sequencename output [option]", process_seq_stat, STR_SEQSTAT, TRUE},
	{"seqsubsky", 2, "seqsubsky sequencename degree [-prefix=]", process_subsky, STR_SEQSUBSKY, TRUE},
	{"set16bits", 0, "set16bits", process_set_32bits, STR_SET16, TRUE},
	{"set32bits", 0, "set32bits", process_set_32bits, STR_SET32, TRUE},
	{"setcompress", 1, "setcompress 0/1 [-type=] [q] [hscale_factor]", process_set_compress, STR_SETCOMPRESS, TRUE},
#ifdef _OPENMP
	{"setcpu", 1, "setcpu number", process_set_cpu, STR_SETCPU, TRUE},
#endif
	{"setext", 1, "setext extension", process_set_ext, STR_SETEXT, TRUE},
	{"setfindstar", 2, "setfindstar sigma roundness", process_set_findstar, STR_SETFINDSTAR, TRUE},
	{"setmag", 1, "setmag magnitude", process_set_mag, STR_SETMAG, FALSE},
	{"setmagseq", 1, "setmagseq magnitude", process_set_mag_seq, STR_SETMAGSEQ, FALSE},
	{"setmem", 1, "setmem ratio", process_set_mem, STR_SETMEM, TRUE},
	{"setref", 2, "setref sequencename image_number", process_set_ref, STR_SETREF, TRUE},
	{"split", 3, "split R G B", process_split, STR_SPLIT, TRUE},
	{"split_cfa", 0, "split_cfa", process_split_cfa, STR_SPLIT_CFA, TRUE},
	{"stack", 1, "stack sequencename [type] [rejection type] [sigma low] [sigma high] [-nonorm, norm=] [-output_norm] [-out=result_filename] [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weighted]", process_stackone, STR_STACK, TRUE},
	{"stackall", 0, "stackall [type] [rejection type] [sigma low] [sigma high] [-nonorm, norm=] [-output_norm] [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weighted]", process_stackall, STR_STACKALL, TRUE},
	{"stat", 0, "stat", process_stat, STR_STAT, TRUE},
	{"subsky", 1, "subsky degree", process_subsky, STR_SUBSKY, TRUE},

	{"threshlo", 1, "threshlo level", process_threshlo, STR_THRESHLO, TRUE},
	{"threshhi", 1, "threshi level", process_threshhi, STR_THRESHHI, TRUE},
	{"thresh", 2, "thresh lo hi", process_thresh, STR_THRESH, TRUE}, /* threshes hi and lo */

	{"unselect", 2, "unselect from to", process_unselect, STR_UNSELECT, FALSE},
	{"unsetmag", 0, "unsetmag", process_unset_mag, STR_UNSETMAG, FALSE},
	{"unsetmagseq", 0, "unsetmagseq", process_unset_mag_seq, STR_UNSETMAGSEQ, FALSE},
	{"unsharp", 2, "unsharp sigma multi", process_unsharp, STR_UNSHARP, TRUE},
	{"visu", 2, "visu low high", process_visu, STR_VISU, FALSE},

	/* wavelet transform in nbr_plan plans */
	{"wavelet", 1, "wavelet nbr_plan type", process_wavelet, STR_WAVELET, TRUE},
	/* reconstruct from wavelet transform and weighs plans with c1, c2, c3... */
	{"wrecons", 2, "wrecons c1 c2 c3 ...", process_wrecons, STR_WRECONS, TRUE},

	{"",0,"",0, STR_NONE, FALSE}
};

#endif /* SRC_CORE_COMMAND_LIST_H_ */
