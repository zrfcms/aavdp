#ifndef __AAVDP_RDF_H__
#define __AAVDP_RDF_H__
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#define THREE_QUARTER_PIINV 0.238732414637843

class RDF
{
public:
    int    numrbin=0, numij=0;
    double *rhoij=nullptr;
    double *rij=nullptr;
    int    **ij=nullptr;
    double **gij=nullptr;
    double dr=0.0;
    RDF(const char *model_path, double rmax, int nbin, bool is_partial=false);
    ~RDF();
    void rdf(const char *rdf_path);
    // void nml(const char *ssf_path);
private:
    double vol;
    void   set_volume(double mat[3][3]);
    void   compute(QB_tools *QB, int natom, double rmax, double rbin, int pairij_id=0);
    void   compute(QB_tools *QB, int natom, double rmax, double rbin, int typei, int typej, int pairij_id);
};

#endif