#ifndef __AAVDP_SND_H__
#define __AAVDP_SND_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#define ZERO_LIMIT 0.000001

struct SND_KNODE{
    int    k[3];
    double theta;
    double intensity;
    int    multiplicity=0;
    SND_KNODE *next=nullptr;
};

//x-ray diffraction
class SND
{
public:
    int    numk=0;
    SND_KNODE  *khead=nullptr;
    SND_KNODE  *ktail=nullptr;
    SND(NMODEL *model, double min2Theta, double max2Theta, double spacing[3], bool is_spacing_auto, bool is_lorentz);
    ~SND();
    void   snd(const char *snd_path, int nbin=0);
private:
    bool   is_lorentz_flag=true;
    double spacingK[3]={0.1, 0.1, 0.1};
    int    kmin[3]={10000, 10000, 10000}, kmax[3]={-10000, -10000, -10000};
    double minTheta=0.0, maxTheta=PI/2.0;
    double intensity_min=1.0e8, intensity_max=0.0;
    void   copy_knode_data(SND_KNODE *knode1, SND_KNODE *knode2);
    void   swap_knode_data(SND_KNODE *knode1, SND_KNODE *knode2);
    void   quick_sort(SND_KNODE *kstart, SND_KNODE *kend);
    void   quick_unique();
    void   add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity);
    void   compute_diffraction_intensity(NMODEL *model);
};

#endif