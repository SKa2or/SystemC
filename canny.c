#include "canny.h"

/*******************************************************************************
* PROCEDURE: canny
* PURPOSE: To perform canny edge detection.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void canny(unsigned char *image, unsigned char *edge, char *fname,
         unsigned char *nms, short int *magnitude,
         short int *delta_x, short int *delta_y, int *tempim,
         short int *smoothedim, int *kernel, int windowsize)
{
    int rows=ROWS;
    int cols=COLS;
    FILE *fpdir=NULL;          /* File to write the gradient image to.     */
    int r, c, pos;

    /****************************************************************************
    * Perform gaussian smoothing on the image using the input standard
    * deviation.
    ****************************************************************************/
    if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
    gaussian_smooth(image, smoothedim, tempim, kernel, windowsize);

    /****************************************************************************
    * Compute the first derivative in the x and y directions.
    ****************************************************************************/
    if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
    derrivative_x_y(smoothedim, delta_x, delta_y);

    /****************************************************************************
    * Compute the magnitude of the gradient.
    ****************************************************************************/
    if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
    magnitude_x_y(delta_x, delta_y, magnitude);
    
    /****************************************************************************
   * Perform non-maximal suppression.
   ****************************************************************************/
   	//if(VERBOSE) printf("Doing the non-maximal suppression.\n");
   	//if((nms = (unsigned char *) calloc(rows*cols,sizeof(unsigned char)))==NULL){
      	//fprintf(stderr, "Error allocating the nms image.\n");
      	//exit(1);
   	//}
    non_max_supp(magnitude, delta_x, delta_y, nms);
    
    /****************************************************************************
   * Use hysteresis to mark the edge pixels.
   ****************************************************************************/
   	//if(VERBOSE) printf("Doing hysteresis thresholding.\n");
   	//if((*edge=(unsigned char *)calloc(rows*cols,sizeof(unsigned char))) ==NULL){
      	//fprintf(stderr, "Error allocating the edge image.\n");
      	//exit(1);
   	//}
    apply_hysteresis(magnitude, nms, edge);


}

/*******************************************************************************
* PROCEDURE: magnitude_x_y
* PURPOSE: Compute the magnitude of the gradient. This is the square root of
* the sum of the squared derivative values.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void magnitude_x_y(short int *delta_x, short int *delta_y,
        short int *magnitude)
{
#pragma HLS INTERFACE ap_bus port = magnitude
    int rows=ROWS;
    int cols=COLS;
    int r, c, pos, sq1, sq2;



    for(r=0,pos=0;r<rows;r++){
        for(c=0;c<cols;c++,pos++){
#pragma HLS PIPELINE II=1
            sq1 = (int)delta_x[pos] * (int)delta_x[pos];
            sq2 = (int)delta_y[pos] * (int)delta_y[pos];
            magnitude[pos] = (sqrt((sq1 + sq2))/1);
            magnitude[pos] = (sqrt((sq1 + sq2)*100)/10);
        }
    }
}

/*******************************************************************************
* PROCEDURE: derrivative_x_y
* PURPOSE: Compute the first derivative of the image in both the x any y
* directions. The differential filters that are used are:
*
*                                          -1
*         dx =  -1 0 +1     and       dy =  0
*                                          +1
*
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void derrivative_x_y(short int *smoothedim,
        short int *delta_x, short int *delta_y)
{
#pragma HLS INTERFACE ap_bus port = smoothedim
#pragma HLS INTERFACE ap_bus port = delta_x
#pragma HLS INTERFACE ap_bus port = delta_y
    int rows=ROWS;
    int cols=COLS;
    int r, c, pos;

    /****************************************************************************
    * Compute the x-derivative. Adjust the derivative at the borders to avoid
    * losing pixels.
    ****************************************************************************/
    if(VERBOSE) printf("   Computing the X-direction derivative.\n");
    for(r=0;r<rows;r++)
    {
        pos = r * cols;
        delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
        pos++;
        for(c=1;c<(cols-1);c++,pos++)
        {
#pragma HLS PIPELINE II=2
            delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
        }
        delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
    }

    /****************************************************************************
    * Compute the y-derivative. Adjust the derivative at the borders to avoid
    * losing pixels.
    ****************************************************************************/
    if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
    for(c=0;c<cols;c++)
    {
        pos = c;
        delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
        pos += cols;
        for(r=1;r<(rows-1);r++,pos+=cols)
        {
#pragma HLS PIPELINE II=2
            delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
        }
        delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
    }
}

/*******************************************************************************
* PROCEDURE: gaussian_smooth
* PURPOSE: Blur an image with a gaussian filter.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void gaussian_smooth(unsigned char *image, short int *smoothedim, int *tempim, int *kernel, int windowsize)
{
    windowsize = WINDOWSIZE;
    int rows=ROWS;
    int cols=COLS;
    int r, c, rr, cc,     /* Counter variables. */
      center;            /* Half of the windowsize. */
    int dot,            /* Dot product summing variable. */
         sum;            /* Sum of the kernel weights variable. */

    /****************************************************************************
    * Create a 1-dimensional gaussian smoothing kernel.
    ****************************************************************************/
    if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
    center = windowsize / 2;

    blur_x(image, center, tempim, kernel);
    blur_y(image, center, tempim, kernel, smoothedim);
}

void blur_x(unsigned char *image, int center, int *tempim, int *kernel)
{
#pragma HLS INTERFACE ap_bus port = image
#pragma HLS INTERFACE ap_bus port = kernel
    center = WINDOWSIZE / 2;
    int rows=ROWS;
    int cols=COLS;
    int r, c, cc;         /* Counter variables. */
    int dot,            /* Dot product summing variable. */
         sum;            /* Sum of the kernel weights variable. */
    /****************************************************************************
    * Blur in the x - direction.
    ****************************************************************************/
    if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
    for(r=0;r<rows;r++)
    {
        for(c=0;c<cols;c++)
        {
#pragma HLS PIPELINE II=5
            dot = 0.0;
            sum = 0.0;
            for(cc=(-center);cc<=center;cc++)
            {
                if(((c+cc) >= 0) && ((c+cc) < cols))
                {
                    dot += image[r*cols+(c+cc)] * kernel[center+cc];
                    sum += kernel[center+cc];
                }
            }
            tempim[r*cols+c] = dot/sum;
        }
    }
}

void blur_y(unsigned char *image, int center, int *tempim,
        int *kernel, short int *smoothedim)
{
#pragma HLS INTERFACE ap_bus port = tempim
    center = WINDOWSIZE / 2;
    int rows=ROWS;
    int cols=COLS;
    int r, c, rr;         /* Counter variables. */
    int dot,            /* Dot product summing variable. */
         sum;            /* Sum of the kernel weights variable. */
    /****************************************************************************
    * Blur in the y - direction.
    ****************************************************************************/
    if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
    for(c=0;c<cols;c++)
    {
        for(r=0;r<rows;r++)
        {
#pragma HLS PIPELINE II=5
            sum = 0.0;
            dot = 0.0;
            for(rr=(-center);rr<=center;rr++)
            {
                if(((r+rr) >= 0) && ((r+rr) < rows))
                {
                   dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
                   sum += kernel[center+rr];
                }
            }
            smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
        }
    }
}

/*******************************************************************************
*******************************************************************************/


/*******************************************************************************
* PROCEDURE: follow_edges
* PURPOSE: This procedure edges is a recursive routine that traces edgs along
* all paths whose magnitude values remain above some specifyable lower
* threshhold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short highval, short lowval)
{
#pragma HLS INTERFACE ap_bus port = edgemapptr
#pragma HLS INTERFACE ap_bus port = edgemagptr
	int cols=COLS;
	short *tempmagptr;
	unsigned char *tempmapptr;
	int i;
	int x[8] = {1,1,0,-1,-1,-1,0,1},
	   y[8] = {0,1,1,1,0,-1,-1,-1};
#pragma HLS ARRAY_PARTITION variable=x complete dim=1
#pragma HLS ARRAY_PARTITION variable=y complete dim=1

	for(i=0;i<8;i++){
		tempmapptr = edgemapptr - y[i]*cols + x[i];
		tempmagptr = edgemagptr - y[i]*cols + x[i];

		if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval))
		{
			*tempmapptr = (unsigned char) EDGE;
			*tempmagptr = highval;
			//follow_edges(tempmapptr,tempmagptr, lowval, cols);
		}
	}
}

/*******************************************************************************
* PROCEDURE: apply_hysteresis
* PURPOSE: This routine finds edges that are above some high threshhold or
* are connected to a high pixel by a path of pixels greater than a low
* threshold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void apply_hysteresis(short int *mag, unsigned char *nms, unsigned char *edge)
{
#pragma HLS INTERFACE ap_bus port = mag
#pragma HLS INTERFACE ap_bus port = nms
#pragma HLS INTERFACE ap_bus port = edge
    int rows=ROWS;
    int cols=COLS;
    int r, c, pos, numedges, lowcount, highcount, lowthreshold, highthreshold,
        i, hist[32768], rr, cc;
    short int maximum_mag, sumpix;

    /****************************************************************************
    * Initialize the edge map to possible edges everywhere the non-maximal
    * suppression suggested there could be an edge except for the border. At
    * the border we say there can not be an edge because it makes the
    * follow_edges algorithm more efficient to not worry about tracking an
    * edge off the side of the image.
    ****************************************************************************/
    for(r=0,pos=0;r<rows;r++){
        for(c=0;c<cols;c++,pos++){
#pragma HLS PIPELINE II=1
            if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
            else edge[pos] = NOEDGE;
        }
    }

    for(r=0,pos=0;r<rows;r++,pos+=cols){
#pragma HLS PIPELINE II=2
        edge[pos] = NOEDGE;
        edge[pos+cols-1] = NOEDGE;
    }
    pos = (rows-1) * cols;
    for(c=0;c<cols;c++,pos++){
#pragma HLS PIPELINE II=2
        edge[c] = NOEDGE;
        edge[pos] = NOEDGE;
    }

    highthreshold = HIGHTHRESHOLD;
    lowthreshold = LOWTHRESHOLD;

    /****************************************************************************
    * This loop looks for pixels above the highthreshold to locate edges and
    * then calls follow_edges to continue the edge.
    ****************************************************************************/
	int total = rows*cols;
	for(pos=0;pos<total;pos++){
		if((edge[pos] != NOEDGE) && (mag[pos] >= highthreshold)){
			edge[pos] = EDGE;
			follow_edges((edge+pos), (mag+pos), highthreshold, lowthreshold);
		}
	}

    /****************************************************************************
    * Set all the remaining possible edges to non-edges.
    ****************************************************************************/
    for(r=0,pos=0;r<rows;r++){
        for(c=0;c<cols;c++,pos++){
        	if(edge[pos] != EDGE){
        		edge[pos] = NOEDGE;
        	}
        }
    }
}

/*******************************************************************************
* PROCEDURE: non_max_supp
* PURPOSE: This routine applies non-maximal suppression to the magnitude of
* the gradient image.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void non_max_supp(short *mag, short *gradx, short *grady,
    unsigned char *result)
{
#pragma HLS INTERFACE ap_bus port = result
#pragma HLS INTERFACE ap_bus port = mag
#pragma HLS INTERFACE ap_bus port = gradx
#pragma HLS INTERFACE ap_bus port = grady
    int nrows=ROWS;
    int ncols=COLS;
    int rowcount, colcount,count;
    short *magrowptr,*magptr;
    short *gxrowptr,*gxptr;
    short *gyrowptr,*gyptr,z1,z2;
    short m00,gx,gy;
    int mag1,mag2,xperp,yperp;
    unsigned char *resultrowptr, *resultptr;


   /****************************************************************************
   * Zero the edges of the result image.
   ****************************************************************************/
    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
        count<ncols; resultptr++,resultrowptr++,count++){
        *resultrowptr = *resultptr = (unsigned char) 0;
    }

    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
        *resultptr = *resultrowptr = (unsigned char) 0;
    }

   /****************************************************************************
   * Suppress non-maximum points.
   ****************************************************************************/
   for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
      gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
      rowcount<=nrows-2;	/* bug fix 3/29/17, RD */
      rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
      resultrowptr+=ncols){
      for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
         resultptr=resultrowptr;colcount<=ncols-2;	/* bug fix 3/29/17, RD */
         colcount++,magptr++,gxptr++,gyptr++,resultptr++){
         m00 = *magptr;
         if(m00 == 0){
            *resultptr = (unsigned char) NOEDGE;
         }
         else{
            xperp = -(gx = *gxptr)*1000/(m00);
            yperp = (gy = *gyptr)*1000/(m00);
         }

         if(gx >= 0){
            if(gy >= 0){
                    if (gx >= gy)
                    {
                        /* 111 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 110 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (gx >= -gy)
                    {
                        /* 101 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 100 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                    }
                }
            }
            else
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {
                        /* 011 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 010 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (-gx > -gy)
                    {
                        /* 001 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 000 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            }

            /* Now determine if the current point is a maximum point */

            if ((mag1 > 1000) || (mag2 > 1000))
            {
                *resultptr = (unsigned char) NOEDGE;
            }
            else
            {
                if (mag2 == 0)
                    *resultptr = (unsigned char) NOEDGE;
                else
                    *resultptr = (unsigned char) POSSIBLE_EDGE;
            }
        }
    }
}
