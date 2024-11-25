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

//electron diffraction
class KED
{
public:
    int    numk=0;
    KED_KNODE  *khead=nullptr;
    KED_KNODE  *ktail=nullptr;
    double Kmagnitude_max=1.70;
    double intensity_min=1.0e8, intensity_max=0.0;
    KED_KNODE *knearest_1=nullptr, *knearest_2=nullptr;
    KED(MODEL *model, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto);
    KED(MODEL *model, int zone[3], double thickness, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto);
    ~KED();
    void   rotate(int x[3], int y[3]);
    void   ked(char *ked_path);
    // void   vtk(char *vtk_path);
    // void   restart(const char *restart_path);
private:
    double axes[3][3]={0.0};
    void   add_k_node(double hkl[3], double K[3], double Kmagnitude, double intensity);
    void   compute_diffraction_intensity(MODEL *model, double threshold, double spacing[3], bool is_spacing_auto);
    void   compute_diffraction_intensity(MODEL *model, int zone[3], double thickness, double spacing[3], bool is_spacing_auto);
    void   filter_diffraction_intensity(double threshold);
    void   find_first_and_second_knearests();
    void   rotate_by_first_knearest(int zone[3]);
    void   img(char *png_path, double *x, double *y, double *value, int num, double limit);
};

#endif