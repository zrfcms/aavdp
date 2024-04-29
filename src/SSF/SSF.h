#ifndef __AAVDP_SSF_H__
#define __AAVDP_SSF_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>
#include "../MATH/MATH.h"

class SSF
{
public:
    int    numqbin=0;
    double *qij=nullptr;
    double *Sij=nullptr;
    SSF(const char *rdf_path, double qmax, double nbin);
    ~SSF();
    void ssf(const char *ssf_path);
private:
    int    numrbin=0;
    double rho0;
    double dr;
    double *rij=nullptr;
    double *gij=nullptr;
    void read(const char *nml_path);
};

#endif