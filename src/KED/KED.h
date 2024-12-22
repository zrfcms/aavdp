#ifndef __AAVDP_KED_H__
#define __AAVDP_KED_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
using namespace std;

#define KED_INTENSITY_LIMIT 1.0e-6
#define KED_KMAG_LIMIT 1.0e-3

#define MAX_LENTH_IN_NAME 50
#define MAX_WORD_NUMBER 100


struct KED_KNODE{
    double hkl[3];
    double K[3];
    double Kmagnitude;
    double intensity;
    KED_KNODE *next=nullptr;
};

extern void copy_knode_data(KED_KNODE *knode1, KED_KNODE *knode2);
extern void swap_knode_data(KED_KNODE *knode1, KED_KNODE *knode2);
extern void quick_sort(KED_KNODE *kstart, KED_KNODE *kend);
extern void quick_sort(KED_KNODE *kstart);

//electron diffraction
class KED
{
public:
    int    numk=0;
    double lambda;
    KED_KNODE  *khead=nullptr;
    KED_KNODE  *ktail=nullptr;
    double Kmagnitude_max=0.0;
    double intensity_min=1.0e8, intensity_max=0.0;
    KED_KNODE *knearest_1=nullptr, *knearest_2=nullptr;
    KED(MODEL *model, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto);
    KED(MODEL *model, int zone[3], double thickness, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto);
    KED(char *ked3_path);
    ~KED();
    void   filter_diffraction_intensity(double threshold);
    void   rotate(double x[3], double y[3]);
    void   ked(char *ked_path);
    void   ked(char *ked_path, double sigma, double dx);
    void   ked3(char *ked3_path);
private:
    double axes[3][3]={0.0};
    void   add_k_node(double hkl[3], double K[3], double Kmagnitude, double intensity);
    void   compute_diffraction_intensity(MODEL *model, double spacing[3], bool is_spacing_auto);
    void   compute_diffraction_intensity(MODEL *model, int zone[3], double thickness, double spacing[3], bool is_spacing_auto);
    void   find_first_and_second_knearests();
    void   rotate_by_first_knearest(int zone[3]);
    void   img(char *png_path, double *x, double *y, double *value, int num, double limit);
};

struct KKD_KNODE{
    double hkl[3];
    double K[3];
    double Kwidth;
    double intensity1;
    double intensity2;
    KKD_KNODE *next=nullptr;
};

//kinematic Kikuchi diffraction
class KKD
{
public:
    KKD(char *ked3_path, double xaxis[3], double yaxis[3], double zaxis[3], double thickness, double threshold, double ratiox, double ratioy, int npx, int npy, char *mode);
    ~KKD();
    int    numk=0;//g
    KKD_KNODE  *khead=nullptr;
    KKD_KNODE  *ktail=nullptr;
    int    numpx=0, numpy=0;//ncol, nrow
    double ***screenK0=nullptr;
    double **screenI=nullptr;
    double thetax=0.0, thetay=0.0;
    double intensity_min=1.0e8, intensity_max=0.0;
    void   kkd(char* kkd_path, char background='b');
    void   kkd(char* kkd_path, double vmax, double vmin, char background='b');
private:
    double axes[3][3]={0.0};
    void   add_k_node(double hkl[3], double K[3], double Kwidth, double intensity1, double intensity2);
    void   compute_Kikuchi_sphere_projection(double xaxis[3], double yaxis[3], double zaxis[3], double kn, double ratiox, double ratioy, int npx, int npy, char *mode);
    void   compute_Kikuchi_intensity_projection(KED *ked, double thickness, double kn);
};

#endif