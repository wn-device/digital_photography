#include <iostream>
#include <tiffio.h>
#include <cstring>
#include <cstdlib>
#include "malloc.h"
#include "fread_dump.h"
#include "hilite_recon.h"
int main (int argc, char **argv)
{
    setlocale (LC_ALL, "");
    setlocale (LC_NUMERIC, "C"); // to set decimal point to "."
    FILE *fp;
    struct imagefloat data_in;
    struct imagefloat * data_in_prt = &data_in;
    const char*imread = "D:/work/data/small_hl.dat";
    const char*imwrtie = "D:/work/data/small_hl_out.dat";
    image_fread(imread,&(data_in_prt));
    // image_dump(imwrtie,data_in_prt);
    int img_size = data_in_prt->height*data_in_prt->width;
    float *red = &(data_in_prt->image[0]);
    float *green = &(data_in_prt->image[img_size]);
    float *blue = &(data_in_prt->image[img_size*2]);
    HLRecovery_inpaint((float *)red,(float *)green,(float *)blue,data_in_prt);
   
}
