#ifndef _AAVDP_DED_H_
#define _AAVDP_DED_H_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "../include/lapacke.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../MODEL/MODEL.h"

using namespace std;

struct DED_BETHE{
    double c1=4.0, c2=8.0, c3=50.0;
    double c_sg=1.0;
};

struct DED_KNODE{
    double k[3];
    double kt[3];
    double kn;
    DED_KNODE *next=nullptr;
};

class DED_KVECTOR
{
public:
    DED_KVECTOR(CELL *cell, double kz[3], double g1[3], double precangle);
    ~DED_KVECTOR();
    int    numk=0;
    double **karray=nullptr;
    double *knarray=nullptr;
private:
    DED_KNODE *khead=nullptr;
    DED_KNODE *ktail=nullptr;
    int    ncircle=10, pcircle=360;
    void   add_k_node(double k[3], double kt[3], double kn);
    void   free_k_node();
};

struct DED_GNODE{
    double hkl[3];//Miller indices
    complex<double> Ug, qg;
    double sg, sig;//excitation error and excitation distance
    bool   is_double_diffrac;
    bool   is_weak, is_strong;
    DED_GNODE *next=nullptr;
    DED_GNODE *nextw=nullptr;//weak
    DED_GNODE *nexts=nullptr;//strong
};

class DED_GVECTOR
{
public:
    DED_GVECTOR(CELL *cell, double kz[3], double fn[3], double cutoff, double precangle, double goffset=0.2);
    ~DED_GVECTOR();
    int    numg=0;
    int    nstrong=0, nweak=0;
    DED_GNODE *ghead=nullptr;
    DED_GNODE *heads=nullptr, *headw=nullptr;
    void   filter_g_node(CELL *cell, DED_BETHE *bethe, double k[3], double fn[3]);
    void   unfilter_g_node();
private:
    DED_GNODE *gtail=nullptr;
    void   add_g_node(double hkl[3], complex<double> Ug, complex<double> qg, double sg, double sig, bool is_double_diffrac);
};

struct DED_INODE{
    double hkl[3];
    double pos[2];
    double intensity;
    DED_INODE *next=nullptr;
};

class DED
{
public:
    DED(CELL *cell, int zone[3], int fnorm[3], double thickness, double dmin, double voltage, double precangle, double c1, double c2, double c3, double c_sg);
    ~DED();
    int    numi=0;
    DED_INODE  *ihead=nullptr;
    DED_INODE  *itail=nullptr;
    void   ded(const char *ded_path);
private:
    double intensity_min=1.0e8, intensity_max=0.0;
    void   compute_dynamic_matrix(complex<double> **dmat, CELL *cell, DED_GVECTOR *gvec);
    void   compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double Z, double KN, int NS);
    void   add_intensity_node(double hkl[3], double pos[2], double intensity);
};

#endif