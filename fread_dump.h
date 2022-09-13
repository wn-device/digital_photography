#ifndef fread_dump
#define fread_dump
#include <stdio.h>
struct imagefloat{
    unsigned int width;
    unsigned int height;
    unsigned int channel;
    float *image;
};
void image_fread(const char *fname,imagefloat ** data_in_prt);
void image_dump(const char *fname,imagefloat *data_out_prt);
void array_to_matrix(float ** red, float ** green,float **blue,float *arr, int row, int col);
#endif