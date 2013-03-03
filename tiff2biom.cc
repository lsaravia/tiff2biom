
/*
 * Convert an RGB tif immage to a chlorophyl-a biomass image
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tiffio.h"
//#include "fortify.h"

#define	MAX_CMAP_SIZE	256

#define	streq(a,b)	(strcmp(a,b) == 0)
#define	strneq(a,b,n)	(strncmp(a,b,n) == 0)

#define	COLOR_DEPTH	8
#define	MAX_COLOR	256

#define	B_DEPTH		5		/* # bits/pixel to use */
#define	B_LEN		(1L<<B_DEPTH)

#define	C_DEPTH		2
#define	C_LEN		(1L<<C_DEPTH)	/* # cells/color to use */

#define	COLOR_SHIFT	(COLOR_DEPTH-B_DEPTH)


float corrFactor[3];
int	bytes_per_pixel;
int	num_colors;
TIFF	*in, *out;
uint32	rowsperstrip = (uint32) -1;
uint16	compression = (uint16) -1;
uint16	bitspersample = 1;
uint16	samplesperpixel;
uint32	imagewidth;
uint32	imagelength;
uint16	predictor = 0;

static	void quant(TIFF*, TIFF*);
static	void usage(void);

#define	CopyField(tag, v) \
	if (TIFFGetField(in, tag, &v)) TIFFSetField(out, tag, v)

int main(int argc, char* argv[])
{
	int i, dither = 0;
	uint16 shortv, config, photometric;
	float floatv;
	uint32 longv;
	int c;
	extern int optind;
	extern char* optarg;


	num_colors = MAX_CMAP_SIZE;
	if (argc < 5)
		usage();

   corrFactor[0]=atof(argv[1]);
   corrFactor[1]=atof(argv[2]);
   corrFactor[2]=atof(argv[3]);

	in = TIFFOpen(argv[4], "r");
	if (in == NULL)
		return (-1);
	TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &imagewidth);
	TIFFGetField(in, TIFFTAG_IMAGELENGTH, &imagelength);
	TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	if (bitspersample != 8 && bitspersample != 16) {
		fprintf(stderr, "%s: Image must have at least 8-bits/sample\n",
		    argv[optind]);
		return (-3);
	}
	if (!TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photometric) ||
	    photometric != PHOTOMETRIC_RGB || samplesperpixel < 3) {
		fprintf(stderr, "%s: Image must have RGB data\n", argv[optind]);
		return (-4);
	}
	TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);
	if (config != PLANARCONFIG_CONTIG) {
		fprintf(stderr, "%s: Can only handle contiguous data packing\n",
		    argv[optind]);
		return (-5);
	}

   /*
	 * STEP 6: scan image, match input values to table entries
	 */
	out = TIFFOpen(argv[5], "w");
	if (out == NULL)
		return (-2);

	CopyField(TIFFTAG_SUBFILETYPE, longv);
	CopyField(TIFFTAG_IMAGEWIDTH, longv);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, (short)COLOR_DEPTH);
	if (compression != (uint16)-1) {
		TIFFSetField(out, TIFFTAG_COMPRESSION, compression);
		switch (compression) {
		case COMPRESSION_LZW:
		case COMPRESSION_DEFLATE:
			if (predictor != 0)
				TIFFSetField(out, TIFFTAG_PREDICTOR, predictor);
			break;
		}
	} else
		CopyField(TIFFTAG_COMPRESSION, compression);

	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, (short)PHOTOMETRIC_MINISBLACK);
	CopyField(TIFFTAG_ORIENTATION, shortv);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, (short)1);
	CopyField(TIFFTAG_PLANARCONFIG, shortv);
	CopyField(TIFFTAG_ROWSPERSTRIP, shortv);
//	TIFFSetField(out, TIFFTAG_ROWSPERSTRIP,
//	    TIFFDefaultStripSize(out, rowsperstrip));
	CopyField(TIFFTAG_MINSAMPLEVALUE, shortv);
	CopyField(TIFFTAG_MAXSAMPLEVALUE, shortv);
	CopyField(TIFFTAG_RESOLUTIONUNIT, shortv);
	CopyField(TIFFTAG_XRESOLUTION, floatv);
	CopyField(TIFFTAG_YRESOLUTION, floatv);
	CopyField(TIFFTAG_XPOSITION, floatv);
	CopyField(TIFFTAG_YPOSITION, floatv);


	quant(in, out);
	(void) TIFFClose(out);
	return (0);
}


char* stuff[] = {
"usage: tiff2biom corrR corrG corrB input.tif output",
NULL
};

static void usage(void)
{
	char buf[BUFSIZ];
	int i;

	setbuf(stderr, buf);
	for (i = 0; stuff[i] != NULL; i++)
		fprintf(stderr, "%s\n", stuff[i]);
	exit(-1);
}




/*
 * straight quantization.  Each pixel is mapped to the colors
 * closest to it.  Color values are rounded to the nearest color
 * table entry.
 */
static void quant(TIFF* in, TIFF* out)
{
	unsigned char	*outline, *inputline;
	register unsigned char	*outptr, *inptr;
	register uint32 i, j;
   uint32 newint[3];

//	register int red, green, blue;
   register unsigned int oldval[3];
   int l,q;
   float pe,ord,foldval;


   q = MAX_COLOR-1;
//   oldval=0;

	inputline = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(in));
	outline = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(out));
	for(i = 0; i < imagelength; i++)
   	{
		if (TIFFReadScanline(in, inputline, i, 0) <= 0)
			break;
		inptr = inputline;
		outptr = outline;

		for (j = 0; j < imagewidth; j++)
		{
			for(l=0; l<3; l++)
			{
			oldval[l] = *inptr++;
			}
	      	foldval = oldval[0]*corrFactor[0]+oldval[1]*corrFactor[1]+oldval[2]*corrFactor[2];
//         	foldval = exp((foldval-847.2664)/-152.4); (Original) Aminot sin corregir R+G+B
         	foldval = exp((foldval-658.6866)/-105.71); //Imagenes de 1999?
//			foldval = oldval[2]*corrFactor[2];
//         	foldval = exp((foldval-233.9509)/-33.5024);  // Banda azul tricromatico 1996
			oldval[0]= floor(foldval+0.5);
			*outptr++ =  oldval[0];
		}

		if (TIFFWriteScanline(out, outline, i, 0) < 0)
//		if (TIFFWriteScanline(out, inputline, i, 0) < 0)
			break;
	}
	_TIFFfree(inputline);
	_TIFFfree(outline);
}


