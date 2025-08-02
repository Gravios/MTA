#include <stdio.h>
#include <unistd.h>
#include "floatfann.h"
#include <stdlib.h>
#include <string.h>


int main( int argc, char* argv[])
{
    const char * netfilename = argv[1];
    const char * datfilename = argv[2];
    const char * resfilename = argv[3];

    FILE* fp_data;
    fp_data = fopen (datfilename,"r");

    int i;
    int s;
    fann_type *calc_out;
    fann_type input[16];
    //float ainput[16];
    char* cinput[16];



    struct fann *ann = fann_create_from_file( netfilename );

    char line [400]; 
    char flc  [400];
    

    FILE* fp_out;
    fp_out = fopen ( resfilename , "w");

    while(fgets(line, 400, fp_data) != NULL)
      {
        line[strcspn(line, "\r\n")] = 0;
        snprintf(flc,sizeof(flc),"%s",line);

	// Read in line of input values
         cinput[0] = strtok(flc," "); 
         input[0] = atof(cinput[0]);
    	 for (i=1; i<16; ++i) { 
    	   cinput[i] = strtok(NULL," "); 
           input[i] = atof(cinput[i]);
    	 } 

	 // run input through network
    	calc_out = fann_run(ann, input); 

        // Write to output file
    	for (i=0; i<6; i++) {
	  fprintf(fp_out,"%f ",calc_out[i]);
	  }
	fprintf(fp_out,"\n"); 

      }

    fclose(fp_data);
    fclose(fp_out);

    fann_destroy(ann);

    return 0;
}
