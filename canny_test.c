#include "canny.h"

int main(int argc, char *argv[])
{
    char *infilename = NULL;  /* Name of the input image */
    char *dirfilename = NULL; /* Name of the output gradient direction image */
    char outfilename[128];    /* Name of the output "edge" image */
    char composedfname[128];  /* Name of the output "direction" image */
    unsigned char *image;     /* The input image */
    unsigned char *edge;      /* The output edge image */
    int rows, cols, i, n, *tempim;           
    float sigma, tlow, thigh;
    unsigned char *nms;
    short int *magnitude, *delta_x, *delta_y; 

    if(argc == 6) dirfilename = infilename;
    else dirfilename = NULL;

    /****************************************************************************
    * Read in the image. This read function allocates memory for the image.
    ****************************************************************************/
    
    for(i=0; i<IMG_NUM; i++)
        {
            n = i % AVAIL_IMG;
            sprintf(infilename, IMG_IN, n+1);
    
			if(VERBOSE) printf("Reading the image %s.\n", infilename);
    		if(read_pgm_image(infilename, &image, ROWS, COLS) == 0)
    		{
        		fprintf(stderr, "Error reading the input image, %s.\n", infilename);
        		exit(1);
   			}

    		rows=ROWS;
    		cols=COLS;

    /****************************************************************************
    * Perform the edge detection. All of the work takes place here.
    ****************************************************************************/
    		if(VERBOSE) printf("Starting Canny edge detection.\n");
    		if(dirfilename != NULL){
        		sprintf(composedfname, "%s_s_%3.2f_l_%3.2f_h_%3.2f.fim", infilename,
        		SIGMA, TLOW, THIGH);
        		dirfilename = composedfname;
    		}
 		
   
    		canny(image, edge, dirfilename, nms, magnitude, delta_x, delta_y, tempim, smoothedim, kernel, windowsize);


    /****************************************************************************
    * Free all of the memory that we allocated except for the edge image that
    * is still being used to store out result.
    ****************************************************************************/
    		//free(smoothedim);
    		free(delta_x);
    		free(delta_y);
    		free(magnitude);
    		free(nms);
    		free(tempim);
    		//free(kernel);




    /****************************************************************************
    * Write out the edge image to a file.
    ****************************************************************************/
    		sprintf(outfilename, "%s_s_%3.2f_l_%3.2f_h_%3.2f.pgm", infilename,
      		SIGMA, TLOW, THIGH);
    		if(VERBOSE) printf("Writing the edge iname in the file %s.\n", outfilename);
    		if(write_pgm_image(outfilename, edge, ROWS, COLS, "", 255) == 0)
    		{
        		fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
        		exit(1);
    		}
    		return(1); /* exit cleanly */
		}
}


/******************************************************************************
* Function: read_pgm_image
* Purpose: This function reads in an image in PGM format. The image can be
* read in from either a file or from standard input. The image is only read
* from standard input when infilename = NULL. Because the PGM format includes
* the number of columns and the number of rows in the image, these are read
* from the file. Memory to store the image is allocated in this function.
* All comments in the header are discarded in the process of reading the
* image. Upon failure, this function returns 0, upon sucess it returns 1.
******************************************************************************/
int read_pgm_image(char *infilename, unsigned char **image, int *rows,
    int *cols)
{
    FILE *fp;
    char buf[71];

    /***************************************************************************
    * Open the input image file for reading if a filename was given. If no
    * filename was provided, set fp to read from standard input.
    ***************************************************************************/
    if(infilename == NULL) fp = stdin;
    else{
        if((fp = fopen(infilename, "r")) == NULL){
            fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
            infilename);
            return(0);
        }
    }

    /***************************************************************************
    * Verify that the image is in PGM format, read in the number of columns
    * and rows in the image and scan past all of the header information.
    ***************************************************************************/
    fgets(buf, 70, fp);
    if(strncmp(buf,"P5",2) != 0){
        fprintf(stderr, "The file %s is not in PGM format in ", infilename);
        fprintf(stderr, "read_pgm_image().\n");
        if(fp != stdin) fclose(fp);
        return(0);
    }
    do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
    sscanf(buf, "%d %d", cols, rows);
    do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

    /***************************************************************************
    * Allocate memory to store the image then read the image from the file.
    ***************************************************************************/
    if(((*image) = (unsigned char *) malloc((*rows)*(*cols))) == NULL){
        fprintf(stderr, "Memory allocation failure in read_pgm_image().\n");
        if(fp != stdin) fclose(fp);
        return(0);
    }
    if((*rows) != fread((*image), (*cols), (*rows), fp)){
        fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
        if(fp != stdin) fclose(fp);
        free((*image));
        return(0);
    }

    if(fp != stdin) fclose(fp);
    return(1);
}




/******************************************************************************
* Function: write_pgm_image
* Purpose: This function writes an image in PGM format. The file is either
* written to the file specified by outfilename or to standard output if
* outfilename = NULL. A comment can be written to the header if coment != NULL.
******************************************************************************/
int write_pgm_image(char *outfilename, unsigned char *image, int rows,
    int cols, char *comment, int maxval)
{
   FILE *fp;

   /***************************************************************************
   * Open the output image file for writing if a filename was given. If no
   * filename was provided, set fp to write to standard output.
   ***************************************************************************/
   if(outfilename == NULL) fp = stdout;
   else{
        if((fp = fopen(outfilename, "w")) == NULL){
            fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
            outfilename);
            return(0);
        }
   }

   /***************************************************************************
   * Write the header information to the PGM file.
   ***************************************************************************/
   fprintf(fp, "P5\n%d %d\n", cols, rows);
   if(comment != NULL)
        if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   /***************************************************************************
   * Write the image data to the file.
   ***************************************************************************/
   if(rows != fwrite(image, cols, rows, fp)){
        fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
        if(fp != stdout) fclose(fp);
        return(0);
   }

   if(fp != stdout) fclose(fp);
   return(1);
}
