#ifndef __AAVDP_SSF_H__
#define __AAVDP_SSF_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <fstream>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../RDF/RDF.h"

class SSF
{
public:
    int    numqbin=0, numij=0;
    double qbin=0.0;
    double *qij=nullptr;
    int    **ij=nullptr;
    double *Smax=nullptr, *Smin=nullptr;
    double **Sij=nullptr;
    SSF(const char *model_path, double qmax, int nbin, bool is_partial=false);
    SSF(RDF *rdf, double qmax, int nbin, bool is_partial=false);
    ~SSF();
    void    ssf(const char *ssf_path, const char *png_path);
private:
    void   compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int pairij_id=0);
    void   compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int typei, int typej, int pairij_id);
};

#endif