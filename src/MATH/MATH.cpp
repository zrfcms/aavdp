#include "MATH.h"

void pseudo_Voigt(double *y, double *x, int num, double x0, double y0, double eta, double w)
{
    if(eta<0.0||eta>1.0){
        printf("[ERROR] Unrecognized mixing parameter for pseudo-Voigt formula\n");
        exit(1);
    }
    double c0=4, c1=4.0*log(2.0);
    double constl=eta*sqrt(c0)/(PI*w), constg=(1.0-eta)*sqrt(c1)/(sqrt(PI)*w);
    for(int i=0;i<num;i++){
        double xk=2.0*(x[i]-x0)/w;
        y[i]=constl/(1.0+c0*xk*xk)+constg*exp(-c1*xk*xk);
        y[i]*=y0;
    }
}

void gaussian(double **value, double *x, double *y, int num, double v0, double x0, double y0, double sigma)
{
    double c0=1.0/(2*PI*sigma*sigma), c1=1.0/(2*sigma*sigma);
    for(int i=0;i<num;i++){
        for(int j=0;j<num;j++){
            value[i][j]=c0*exp(-c1*((y[i]-y0)*(y[i]-y0)+(x[j]-x0)*(x[j]-x0)));
            value[i][j]*=v0;
        }
    }
}

void vector_rotate(double r_v[3], int R[3][3], double v[3])
{
    double c_v[3];
    vector_copy(c_v, v);
    for(int i=0;i<3;i++){
        r_v[i]=0.0;
        for(int j=0;j<3;j++){
            r_v[i]+=(double)R[i][j]*c_v[j];
        }
    }
}