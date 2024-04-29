#ifndef __AAVDP_RDF_H__
#define __AAVDP_RDF_H__
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"

class RDF
{
public:
    int    numij=0, numrbin=0;
    int    **ij=nullptr;
    double *rho0ij=nullptr;
    double *rij=nullptr;
    double **gij=nullptr;
    RDF(const char *model_path, double rmax, int nbin, bool is_partial_flag);
    ~RDF();
    void rdf(const char *rdf_path);
    void nml(const char *ssf_path);
private:
    double dr;
    double vol;
    void   volume(QB_tools *QB);
    void   compute(QB_tools *QB, double rmax);
    void   compute(QB_tools *QB, double rmax, int typei, int typej, int pair_id);
};

#endif