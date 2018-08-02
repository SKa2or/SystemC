//Final
//Shuo Jiang

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define VERBOSE 0

//#include "systemc.h"

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define ROWS 1520
#define COLS 2704
#define SIZE COLS*ROWS
#define WINDOWSIZE 21

#define VIDEONAME "Engineering"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */

#define MAX_INT ((unsigned)(-1)>>1)
#define MIN_INT (~MAX_INT)

#define MAX_SHORT_INT ((unsigned short)(-1)>>1)
#define MIN_SHORT_INT (~MAX_SHORT_INT)

#define HIGHTHRESHOLD 2803
#define LOWTHRESHOLD 841


int read_pgm_image(char *infilename, unsigned char **image, int *rows, int *cols);
int write_pgm_image(char *outfilename, unsigned char *image, int rows, int cols, char *comment, int maxval);

void canny(unsigned char *image, unsigned char *edge, char *fname, unsigned char *nms, short int *magnitude,
    short int *delta_x, short int *delta_y, int *tempim,short int *smoothedim, int *kernel, int windowsize);
void gaussian_smooth(unsigned char *image, short int *smoothedim, int *tempim, int *kernel, int windowsize);
void derrivative_x_y(short int *smoothedim, short int *delta_x, short int *delta_y);
void magnitude_x_y(short int *delta_x, short int *delta_y, short int *magnitude);
void apply_hysteresis(short int *mag, unsigned char *nms, unsigned char *edge);
double angle_radians(double x, double y);
void non_max_supp(short *mag, short *gradx, short *grady, unsigned char *result);
void blur_x(unsigned char *image, int center, int *tempim, int *kernel);
void blur_y(unsigned char *image, int center, int *tempim, int *kernel, short int *smoothedim);
void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short highval, short lowval);

/*typedef struct Image_s
{
    unsigned char img[SIZE];
    
    Image_s(void)
    {
        for (int i=0; i<SIZE; i++)
        {
            img[i] = 0;
        }
    }
    
    Image_s& operator=(const Image_s& copy)
    {
        for (int i=0; i<SIZE; i++)
        {
            img[i] = copy.img[i];
        }
        return *this;
    }
    
    operator unsigned char*()
    {
        return img;
    }
    
    unsigned char& operator[](const int index)
    {
        return img[index];
    }
} IMAGE;

typedef struct Simage_s
{
    short int simg[SIZE];
    
    Simage_s(void)
    {
        for (int i=0; i<SIZE; i++)
        {
            simg[i] = 0;
        }
    }
    
    Simage_s& operator=(const Simage_s& copy)
    {
        for (int i=0; i<SIZE; i++)
        {
            simg[i] = copy.simg[i];
        }
        return *this;
    }
    
    operator short int*()
    {
        return simg;
    }
    
    short int& operator[](const int index)
    {
        return simg[index];
    }
} SIMAGE;

typedef struct Fkernal_s
{
    float fkernal[WINSIZE];
    
    Fkernal_s(void)
    {
        for (int i=0; i<WINSIZE; i++)
        {
            fkernal[i] = 0;
        }
    }
    
    Fkernal_s& operator=(const Fkernal_s& copy)
    {
        for (int i=0; i<WINSIZE; i++)
        {
            fkernal[i] = copy.fkernal[i];
        }
        return *this;
    }
    
    operator float*()
    {
        return fkernal;
    }
    
    float& operator[](const int index)
    {
        return fkernal[index];
    }
} FKERNAL;

typedef struct Fimage_s
{
    float fimg[SIZE];
    
    Fimage_s(void)
    {
        for (int i=0; i<SIZE; i++)
        {
            fimg[i] = 0;
        }
    }
    
    Fimage_s& operator=(const Fimage_s& copy)
    {
        for (int i=0; i<SIZE; i++)
        {
            fimg[i] = copy.fimg[i];
        }
        return *this;
    }
    
    operator float*()
    {
        return fimg;
    }
    
    float& operator[](const int index)
    {
        return fimg[index];
    }
} FIMAGE;*/


