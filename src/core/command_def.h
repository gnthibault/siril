#ifndef SRC_CORE_COMMAND_DEF_H_
#define SRC_CORE_COMMAND_DEF_H_

#define STR_NONE ""

#define STR_ADDMAX N_("Computes a new image IMG with IMG_1 and IMG_2. The pixel of IMG_1 is replaced by the pixel at the same coordinates of IMG_2 if the intensity of 2 is greater than 1")
#define STR_ASINH N_("ASINH command stretches the image for show faint objects, while simultaneously, preserve the structure of bright objects of the field")

#define STR_BG N_("Returns the background level of the image loaded in memory")
#define STR_BGNOISE N_("Returns the background noise level")

#define STR_CD N_("Set the new current working directory. The argument \"directory\" can contain the ~ token, expanded as the home directory, directories with spaces in the name can be protected using single or double quotes")
#define STR_CDG N_("Returns the coordinates of the center of gravity of the image")
#define STR_CLAHE N_("Equalizes the histogram of an image using Contrast Limited Adaptive Histogram Equalization")
#define STR_CLEAR N_("Clears the graphical output logs")
#define STR_CLEARSTAR N_("Clear all the stars saved in memory and displayed on the screen")
#define STR_CLOSE N_("Properly closes the opened image and the opened sequence, if any")
#define STR_CONVERT N_("Convert all images in a known format into Siril's FITS images. The argument \"basename\" is the basename of the new sequence. For FITS images, Siril will try to make a symbolic link. If not possible, files are copied.  The option \"-debayer\" applies demosaicing to images, in this case no symbolic link are done. \"-start=index\" sets the starting index parameter and the \"-out=\" option converts files into another directory" )
#define STR_CONVERTRAW N_("Convert DSLR RAW files into Siril's FITS images. The argument \"basename\" is the basename of the new sequence. The option \"-debayer\" applies demosaicing to images while \"-fitseq\" is an option that sets the output file type to a FITS sequence. \"-start=index\" sets the starting index parameter. It is possible to convert files in another directory with the \"-out=\" option")
#define STR_COSME N_("Apply the local mean to a set of pixels on the in-memory image (cosmetic correction). The coordinates of these pixels are in an ASCII file [.lst file]. COSME is adapted to correct residual hot and cold pixels after preprocessing")
#define STR_COSME_CFA N_("Same function as COSME but applying to RAW CFA images")
#define STR_CROP N_("It can be used with the GUI: if a selection has been made with the mouse, calling the CROP command without arguments crops it on this selection. Otherwise, or in scripts, arguments have to be given, with \"x\" and \"y\" being the coordinates of the top left corner, and \"width\" and \"height\" the size of the selection")

#define STR_DDP N_("Performs a DDP (digital development processing) as described first by Kunihiko Okano. This implementation is the one described in IRIS. It combines a linear distribution on low levels (below \"level\") and a non-linear on high levels. It uses a Gaussian filter of sigma \"sigma\" and multiplies the resulting image by \"coef\". The typical values for \"sigma\" are included between 0.7 and 2")

#define STR_ENTROPY N_("Computes the entropy of the opened image on the displayed layer, only in the selected area if one has been selected or in the whole image. The entropy is one way of measuring the noise or the details in an image")
#define STR_EXIT N_("Quits the application")
#define STR_EXTRACT N_("Extracts \"NbPlans\" planes of wavelet domain")
#define STR_EXTRACTHA N_("Extracts Ha signal from a CFA image. The output file name starts with the prefix \"Ha_\"")
#define STR_EXTRACTHAOIII N_("Extracts Ha and OIII signals from a CFA image. The output file name start with the prefix \"Ha_\" and \"OIII_\"")

#define STR_FDIV N_("Divides the image in memory by the image given in argument. The resulting image is multiplied by the value of the \"scalar\" argument. See also IDIV")
#define STR_FFTD N_("Applies a Fast Fourier Transform to the image loaded in memory. \"Modulus\" and \"phase\" given in argument are saved in FITS files")
#define STR_FFTI N_("Retrieves corrected image applying an inverse transformation. The \"modulus\" and \"phase\" used are the files given in argument")
#define STR_FILL N_("Fills the whole current image (or selection) with pixels having the \"value\" intensity")
#define STR_FILL2 N_("Same command as FILL but this is a symmetric fill of a region defined by the mouse. Used to process an image in the Fourier (FFT) domain")
#define STR_FIND_COSME N_("Applies an automatic detection of cold and hot pixels following the thresholds written in arguments")
#define STR_FIND_COSME_CFA N_("Same command as FIND_COSME but for monochromatic CFA images")
#define STR_FIND_HOT N_("Provides a list file \"filename\" (format text) in the working directory which contains the coordinates of the pixels which have an intensity \"hot_sigma\" times higher and \"cold_sigma\" lower than standard deviation. We generally use this command on a master-dark file")
#define STR_FINDSTAR N_("Detects stars having a level greater than a threshold computed by Siril. The algorithm is based on the publication of Mighell, K. J. 1999, in ASP Conf. Ser., Vol. 172, Astronomical Data Analysis Software and Systems VIII, eds. D. M. Mehringer, R. L. Plante, and D. A. Roberts (San Francisco: ASP), 317. After that, a PSF is applied and Siril rejects all detected structures that don't fulfill a set of prescribed detection criteria. Finally, a circle is drawn around detected stars. See also the command CLEARSTAR")
#define STR_FMEDIAN N_("Performs a median filter of size \"ksize\" x \"ksize\" (\"ksize\" MUST be odd) to the original image with a modulation parameter \"modulation\". The output pixel is computed as : out=mod x m + (1 − mod) x in, where m is the median-filtered pixel value. A modulation's value of 1 will apply no modulation")
#define STR_FMUL N_("Multiplies the loaded image by the \"scalar\" given in argument")
#define STR_FIXBANDING N_("Tries to remove the canon banding. Argument \"amount\" define the amount of correction. \"Sigma\" defines a protection level of the algorithm, higher sigma gives higher protection")

#define STR_GAUSS N_("Performs a Gaussian filter with the given \"sigma\"")
#define STR_GREY_FLAT N_("The function equalizes the mean intensity of RGB layers in a CFA images")

#define STR_HELP N_("Gives the available commands")
#define STR_HISTO N_("Calculates the histogram of the image channel in memory and produces file histo_[channel name].dat in the working directory")

#define STR_IADD N_("Adds the image in memory to the image given in argument")
#define STR_IDIV N_("Divides the image in memory by the image given in argument. See also FDIV")
#define STR_IMUL N_("Multiplies the image in memory by the image given in argument")
#define STR_ISUB N_("Subtracts the image in memory by the image given in argument")

#define STR_LINK N_("Link all FITS images in the working directory with the basename given in argument. If no symbolic link could be created, files are copied. It is possible to convert files in another directory with the \"-out=\" option")
#define STR_LOAD N_("Loads the image \"filename\"; it first attempts to load \"filename\", then \"filename\".fit and finally \"filename\".fits and after, all supported format, aborting if none of these are found. These scheme is applicable to every Siril command implying reading files. Fits headers MIPS-HI and MIPS-LO are read and their values given to the current viewing levels. Writing a known extension at the end of \"filename\" will load the image \"filename\".ext: this is used when numerous files have the same name but not the same extension")
#define STR_LOG N_("Computes and applies a logarithmic scale to the current image")
#define STR_LMATCH N_("Computes a linear function between a reference image and a target. The function is then applied to the target image to match it to the reference one. The algorithm will ignore all reference pixels whose values are outside of the [\"low\", \"high\"] range")
#define STR_LS N_("Lists files and directories in the working directory")

#define STR_MIRRORX N_("Rotates the image around a vertical axis")
#define STR_MIRRORY N_("Rotates the image around an horizontal axis")
#define STR_MTF N_("Applies midtone transfer function to the current loaded image")

#define STR_NEG N_("Shows the negative view of the current image")
#define STR_NEW N_("Creates a new image filled with zeros with a size of \"width\" x \"height\". The image is in 16-bit format, and it contains \"nb_channel\" channels, \"nb_channel\" being 1 or 3. It is not saved, but displayed and can be saved afterwards")
#define STR_NOZERO N_("Replaces null values by \"level\" values. Useful before an idiv or fdiv operation")

#define STR_OFFSET N_("Adds the constant \"value\" to the current image. This constant can take a negative value. As Siril uses unsigned FITS files, if the intensity of the pixel become negative its value is replaced by 0 and by 65535 (for a 16-bit file) if the pixel intensity overflows")

#define STR_PREPROCESS N_("Preprocesses the sequence \"sequencename\" using bias, dark and flat given in argument. It is possible to specify if images are CFA for cosmetic correction purposes with the option \"-cfa\" and also to demosaic images at the end of the process with \"-debayer\". The \"-fix_xtrans\" option is dedicated to X-Trans files by applying a correction on darks and biases to remove an ugly square pattern and the \"-equalize_cfa\" option equalizes the mean intensity of RGB layers of the CFA flat master. It is also possible to optimize the dark subtraction with \"-opt\". The output sequence name starts with the prefix \"pp_\" unless otherwise specified with option \"-prefix=\". If \"-fitseq\" is provided, the output sequence will be a FITS sequence (single file).\n\nNote that only hot pixels are corrected in cosmetic correction process")
#define STR_PSF N_("Performs a PSF (Point Spread Function) on the selected star")

#define STR_REGISTER N_("Performs geometric transforms on images of the sequence given in argument so that they may be superimposed on the reference image. The output sequence name starts with the prefix \"r_\" unless otherwise specified with \"-prefix=\" option. Using stars for registration, this algorithm only works with deepsky images. The registration is done on the green layer for RGB images. The option \"-drizzle\" activates the sub-pixel stacking, either by up-scaling by 2 the images created in the rotated sequence or by setting a flag that will proceed to the up-scaling during stacking if \"-norot\" is passed")
#define STR_RESAMPLE N_("Resamples image with a factor \"factor\"")
#define STR_RELOADSCRIPTS N_("Rescans the scripts folders and updates scripts menu")
#define STR_RGRADIENT N_("Creates two images, with a radial shift (\"dR\" in pixels) and a rotational shift (\"dalpha\" in degrees) with respect to the point (\"xc\", \"yc\"). Between these two images, the shifts have the same amplitude, but an opposite sign. The two images are then added to create the final image. This process is also called Larson Sekanina filter")
#define STR_RL N_("Restores an image using the Richardson-Lucy method. \"Sigma\" is the size of the kernel to be applied, while \"corner_radius_boost\" is a value which is added to Gaussian sigma for the tiles in the corners of an image. \"Iterations\" is the number of iterations to be performed")
#define STR_RMGREEN N_("Chromatic noise reduction filter. It removes green noise in the current image. This filter is based on PixInsight's SCNR Average Neutral algorithm and it is the same filter used by HLVG plugin in Photoshop. \"Type\"=1 stands for Average Neutral Protection, while \"type\"=2 stands for Maximum Neutral Protection")
#define STR_ROTATE N_("Rotates the image of an angle of \"degree\" value")
#define STR_ROTATEPI N_("Rotates the image of an angle of 180° around its center. This is equivalent to the command \"ROTATE 180\" or \"ROTATE -180\"")

#define STR_SATU N_("Enhances the global saturation of the image. Try iteratively to obtain best results")
#define STR_SAVE N_("Saves current image to \"filename\".fit. Fits headers MIPS-HI and MIPS-LO are added with values corresponding to the current viewing levels")
#define STR_SAVEBMP N_("Saves current image under the form of a bitmap file with 8-bit per channel: \"filename\".bmp (BMP 24-bit).")
#define STR_SAVEJPG N_("Saves current image into a JPG file: \"filename\".jpg. You have the possibility to adjust the quality of the compression. A value 100 for \"quality\" parameter offers best fidelity while a low value increases the compression ratio. If no value is specified, it holds a value of 100")
#define STR_SAVEPNG N_("Saves current image into a PNG file: \"filename\".png")
#define STR_SAVEPNM N_("Saves current image under the form of a Netpbm file format with 16-bit per channel. The extension of the output will be \"filename\".ppm for RGB image and \"filename\".pgm for gray-level image")
#define STR_SAVETIF N_("Saves current image under the form of a uncompressed TIFF file with 16-bit per channel: \"filename\".tif")
#define STR_SAVETIF32 N_("Same command as SAVE_TIF but the output file is saved in 32-bit per channel: \"filename\".tif")
#define STR_SAVETIF8 N_("Same command as SAVE_TIF but the output file is saved in 8-bit per channel: \"filename\".tif")
#define STR_SELECT N_("This command allows easy mass selection of images in the loaded sequence (from \"from\" to \"to\" included)")
#define STR_SEQCROP N_("Crops the loaded sequence. The output sequence name starts with the prefix \"cropped_\" unless otherwise specified with \"-prefix=\" option")
#define STR_SEQEXTRACTHA N_("Same command as EXTRACT_HA but for the sequence \"sequencename\". The output sequence name starts with the prefix \"Ha_\" unless otherwise specified with option \"-prefix=\"")
#define STR_SEQEXTRACTHAOIII N_("Same command as EXTRACT_HAOIII but for the sequence \"sequencename\". The output sequence name start with the prefix \"Ha_\" and \"OIII_\"")
#define STR_SEQFIND_COSME N_("Same command as FIND_COSME but for the sequence \"sequencename\". The output sequence name starts with the prefix \"cc_\" unless otherwise specified with \"-prefix=\" option")
#define STR_SEQFIND_COSME_CFA N_("Same command as FIND_COSME_CFA but for the sequence \"sequencename\". The output sequence name starts with the prefix \"cc_\" unless otherwise specified with \"-prefix=\" option")
#define STR_SEQMTF N_("Same command as MTF but for the sequence \"sequencename\".  The output sequence name starts with the prefix \"mtf_\" unless otherwise specified with \"-prefix=\" option")
#define STR_SEQPSF N_("Same command as PSF but works for sequences. Results are dumped in the console in a form that can be used to produce brightness variation curves")
#define STR_SEQSPLIT_CFA N_("Same command as SPLIT_CFA but for the sequence \"sequencename\". The output sequence name starts with the prefix \"CFA_\" unless otherwise specified with \"-prefix=\" option")
#define STR_SEQSUBSKY N_("Same command as SUBSKY but for the sequence \"sequencename\".  The output sequence name starts with the prefix \"bkg_\" unless otherwise specified with \"-prefix=\" option")
#define STR_SET16 N_("Disallow images to be saved with 32 bits per channel on processing, use 16 instead")
#define STR_SET32 N_("Allow images to be saved with 32 bits per channel on processing")
#define STR_SETCOMPRESS N_("Defines if images are compressed or not: 0 means no compression. If compression is enabled, the type must be explicitly written in the option \"-type=\" (\"rice\", \"gzip1\", \"gzip2\" and \"hcompress\"). Associated to the compression, the quantization value must follow [0, 256] and in the case of \"hcompress\", the hcompress scale factor [0, 256] as well. For example, \"setcompress 1 -type=rice 16\" set the rice compression with a quantization of 16")
#define STR_SETCPU N_("Defines the number of processing threads used for calculation. Can be as high as the number of virtual threads existing on the system, which is the number of CPU cores or twice this number if hyperthreading (Intel HT) is available")
#define STR_SETEXT N_("Sets the extension used and recognized by sequences. The argument \"extension\" can be \"fit\", \"fts\" or \"fits\"")
#define STR_SETFINDSTAR N_("Defines thresholds above the noise and star roundness for stars detection with FINDSTAR and REGISTER commands. \"Sigma\" must be greater or equal to 0.05 and \"roundness\" between 0 and 0.9")
#define STR_SETMAG N_("Calibrates the magnitude by selecting a star and giving the known apparent magnitude. All PSF computations will return the calibrated apparent magnitude afterwards, instead of an apparent magnitude relative to ADU values. To reset the magnitude constant see UNSETMAG")
#define STR_SETMAGSEQ N_("This command is only valid after having run SEQPSF or its graphical counterpart (select the area around a star and launch the PSF analysis for the sequence, it will appear in the graphs). This command has the same goal as SETMAG but recomputes the reference magnitude for each image of the sequence where the reference star has been found. When running the command, the last star that has been analysed will be considered as the reference star. Displaying the magnitude plot before typing the command makes it easy to understand. To reset the reference star and magnitude offset, see UNSETMAGSEQ")
#define STR_SETMEM N_("Sets a new ratio of free memory on memory used for stacking. Value should be between 0.05 and 2, depending on other activities of the machine. A higher ratio should allow siril to stack faster, but setting the ratio of memory used for stacking above 1 will require the use of on-disk memory, which is very slow and unrecommended")
#define STR_SETREF N_("Sets the reference image of the loaded sequence")
#define STR_SPLIT N_("Splits the color image into three distinct files (one for each color) and save them in \"r\" \"g\" and \"b\" file")
#define STR_SPLIT_CFA N_("Splits the CFA image into four distinct files (one for each channel) and save them in files")
#define STR_STACK N_("Stacks the \"sequencename\" sequence, using options. The allowed types are: sum, max, min, med or median, and rej or mean that requires the use of additional arguments \"sigma low\" and \"high\" used for the Winsorized sigma clipping rejection algorithm (cannot be changed from here).\nDifferent types of normalisation are allowed: \"-norm=add\" for addition, \"-norm=mul\" for multiplicative. Options \"-norm=addscale\" and \"-norm=mulscale\" apply same normalization but with scale operations. Finally, \"-output_norm\" applies a normalization at the end of the stacking to rescale result in the [0, 1] range.\nIf no argument other than the sequence name is provided, sum stacking is assumed.\nResult image's name can be set with the \"-out=\" option.\nStacked images can be selected based on some filters, like manual selection or best FWHM, with some of the \"-filter-*\" options.\nSee the command reference for the complete documentation on this command")
#define STR_STACKALL N_("Opens all sequences in the CWD and stacks them with the optionally specified stacking type and filtering or with sum stacking. See STACK command for options description")
#define STR_STAT N_("Returns global statistics of the current image. If a selection is made, the command returns statistics within the selection")
#define STR_SUBSKY N_("Computes the level of the local sky background thanks to a polynomial function of an order ''degree'' and subtracts it from the image. A synthetic image is then created and subtracted from the original one")

#define STR_THRESHLO N_("Replaces values below \"level\" with \"level\"")
#define STR_THRESHHI N_("Replaces values above \"level\" with \"level\"")
#define STR_THRESH N_("Replaces values below \"lo\" with \"lo\" and values above \"hi\" with \"hi\"")

#define STR_UNSELECT N_("Allows easy mass unselection of images in the loaded sequence (from \"from\" to \"to\" included). See SELECT")
#define STR_UNSETMAG N_("Reset the magnitude calibration to 0. See SETMAG")
#define STR_UNSETMAGSEQ N_("Resets the magnitude calibration and reference star for the sequence. See SETMAGSEQ")
#define STR_UNSHARP N_("Applies to the working image an unsharp mask with sigma \"sigma\" and coefficient \"multi\"")

#define STR_VISU N_("Displays an image with \"low\" and \"high\" as the low and high threshold")

#define STR_WAVELET N_("Computes the wavelet transform on \"nbr_plan\" plans using linear (type=1) or bspline (type=2) version of the 'à trous' algorithm. The result is stored in a file as a structure containing the planes, ready for weighted reconstruction with WRECONS")
#define STR_WRECONS N_("Reconstructs to current image from the planes previously computed with wavelets and weighted with coefficients \"c1\", \"c2\", ..., \"cn\" according to the number of planes used for wavelet transform")



#endif /* SRC_CORE_COMMAND_DEF_H_ */
