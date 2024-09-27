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
    double c_sg=0.2;
};

struct DED_GNODE{
    double hkl[3];//Miller indices
    complex<double> Ug, qg;
    double sg;//excitation error
    double xg;// extinction distance
    bool   is_double_diffrac;
    bool   is_weak, is_strong;
    DED_GNODE *next=nullptr;
    DED_GNODE *nextw=nullptr;//weak
    DED_GNODE *nexts=nullptr;//strong
};

class DED_GVECTOR
{
public:
    DED_GVECTOR(CELL *cell, DED_BETHE *bethe, double kz[3], double fn[3], double cutoff);
    ~DED_GVECTOR();
    int    numg=0;
    int    nstrong=0, nweak=0;
    DED_GNODE *ghead=nullptr;
    DED_GNODE *heads=nullptr, *headw=nullptr;
private:
    DED_GNODE *gtail=nullptr;
    void   add_g_node(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac);
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
    DED(CELL *cell, DED_BETHE *bethe, int zone[3], int fnorm[3], double voltage, double thickness, double dmin);
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