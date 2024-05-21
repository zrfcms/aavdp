#ifndef __AAVDP_SSF_H__
#define __AAVDP_SSF_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <fstream>
#include <ctime>
#include "../MATH/MATH.h"
#include "../RDF/RDF.h"

class SSF
{
public:
    int    numqbin=0, numij=0;
    double *qij=nullptr;
    int    **ij=nullptr;
    double **Sij=nullptr;
    SSF(const char *model_path, double qmax, int nbin, bool is_partial=false);
    ~SSF();
    void    ssf(const char *ssf_path);
private:
    void   compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int pairij_id=0);
    void   compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int typei, int typej, int pairij_id);
};

#endif