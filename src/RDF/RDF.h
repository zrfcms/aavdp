#ifndef __AAVDP_RDF_H__
#define __AAVDP_RDF_H__
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
#define THREE_QUARTER_PIINV 0.238732414637843

extern void int_to_str(char str[], int num);

class RDF
{
public:
    int    numrbin=0, numij=0;
    int    **ij=nullptr;
    double *rhoij=nullptr;
    double *rij=nullptr;
    double **gij=nullptr;
    double rbin=0.0;
    double *gmax=nullptr, *gmin=nullptr;
    RDF(char *model_path, double rmax, int nbin, bool is_partial=false);
    ~RDF();
    void rdf(char *rdf_path);
private:
    double vol;
    void   set_volume(double mat[3][3]);
    void   compute(QB_tools *QB, int natom, double rmax, int pairij_id=0);
    void   compute(QB_tools *QB, int natom, double rmax, int typei, int typej, int pairij_id);
    void   img(char *png_path, double *x, double *y, int num, double ymin, double ymax);
};

class SSF
{
public:
    int    numqbin=0, numij=0;
    double qbin=0.0;
    double *qij=nullptr;
    int    **ij=nullptr;
    double *Smax=nullptr, *Smin=nullptr;
    double **Sij=nullptr;
    SSF(char *model_path, double qmax, int nbin, bool is_partial=false);
    SSF(RDF *rdf, double qmax, int nbin, bool is_partial=false);
    ~SSF();
    void    ssf(char *ssf_path);
private:
    void   compute(QB_tools *QB, double qmax, double spacingK[3], int NspacingK[3], int pairij_id=0);
    void   compute(QB_tools *QB, double qmax, double spacingK[3], int NspacingK[3], int typei, int typej, int pairij_id);
    void   img(char *png_path, double *x, double *y, int num, double ymin, double ymax);
};

#endif