#ifndef hilite_recon
#define hilite_recon

struct Imatrices{
    double xyz_cam[3][3]={0};
};
void HLRecovery_inpaint(float* red, float* green, float* blue, struct imagefloat *imagedata);
#endif
