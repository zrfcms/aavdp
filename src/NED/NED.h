#ifndef __AAVDP_NED_H__
#define __AAVDP_NED_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"

struct NED_KNODE{
    int    k[3];
    double theta;
    double intensity;
    int    multiplicity=0;
    NED_KNODE *next=nullptr;
};

//neutron diffraction
class NED
{
public:
    int    numk=0;
    NED_KNODE  *khead=nullptr;
    NED_KNODE  *ktail=nullptr;
    NED(NMODEL *model, double min2Theta, double max2Theta, int lp_type, double spacing[3], bool is_spacing_auto);
    NED(char *ned_path);
    ~NED();
    void   ned(char *ned_path, char *png_path);
    void   ned(char *ned_path, char *png_path, double peak_parameter, double scherrer_lambda, double scherrer_size, double dt);
private:
    double spacingK[3]={0.1, 0.1, 0.1};
    int    kmin[3]={10000, 10000, 10000}, kmax[3]={-10000, -10000, -10000};
    double minTheta=0.0, maxTheta=PI/2.0;
    double intensity_min=1.0e8, intensity_max=0.0;
    void   copy_knode_data(NED_KNODE *knode1, NED_KNODE *knode2);
    void   swap_knode_data(NED_KNODE *knode1, NED_KNODE *knode2);
    void   quick_sort(NED_KNODE *kstart, NED_KNODE *kend);
    void   quick_unique();
    void   add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity);
    void   compute_diffraction_intensity(NMODEL *model, int lp_type);
};

#endif