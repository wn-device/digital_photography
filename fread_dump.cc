#include <stdio.h>
#include "malloc.h"
#include "fread_dump.h"

void image_fread(const char *fname,imagefloat ** data_in_prt){
FILE *fp  = fopen(fname,"rb");
fread(&((*data_in_prt)->width),sizeof(unsigned int),1,fp);
fread(&((*data_in_prt)->height),sizeof(unsigned int),1,fp);
fread(&((*data_in_prt)->channel),sizeof(unsigned int),1,fp);
int image_size_byte = (*data_in_prt)->width*(*data_in_prt)->height*(*data_in_prt)->channel*sizeof(float);
(*data_in_prt)->image = (float *)malloc(image_size_byte);
fread((*data_in_prt)->image,image_size_byte,1,fp);
fclose(fp);
printf("fread correct\n");
}
void image_dump(const char *fname,imagefloat *data_out_prt){
FILE *fp    = fopen(fname,"wb");
fwrite(&(data_out_prt->width),sizeof(unsigned int),1,fp);
fwrite(&(data_out_prt->height),sizeof(unsigned int),1,fp);
fwrite(&(data_out_prt->channel),sizeof(unsigned int),1,fp);
int image_size_byte = data_out_prt->width*data_out_prt->height*data_out_prt->channel*sizeof(float);
fwrite(data_out_prt->image,image_size_byte,1,fp);
fclose(fp);
printf("fwrite correct\n");
}