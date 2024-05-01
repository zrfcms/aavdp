#include "XRD.h"

XRD::XRD(MODEL *model, double spacingK[3], double min2Theta, double max2Theta, bool is_lorentz, bool is_spacing_auto)
{
    printf("Starting...\n");
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacingK);
    }
    double minTheta=min2Theta/2.0*DEG_TO_RAD, maxTheta=max2Theta/2.0*DEG_TO_RAD;
    double Kmagnitude_max=2.0*sin(maxTheta)/model->lambda;
    for(int i=0;i<3;i++){
        int NspacingK=ceil(Kmagnitude_max/spacingK[i]);
        kmin[i]=-NspacingK; kmax[i]=NspacingK;
    }
    is_lorentz_flag=is_lorentz;
    count_diffraction_vector(model, spacingK, minTheta, maxTheta);
    compute_diffraction_intensity(model, spacingK, minTheta, maxTheta);
    printf("Ending...\n");
}

XRD::XRD(MODEL *model, bool is_lorentz)
{
    printf("Starting...\n");
    double Kmagnitude_max[3];
    vector_constant(Kmagnitude_max, 2.0/model->lambda, model->dimension);
    for(int i=0;i<3;i++){
        int NspacingK=ceil(Kmagnitude_max[i]);
        kmin[i]=-NspacingK; kmax[i]=NspacingK;
    }
    is_lorentz_flag=is_lorentz;
    count_diffraction_vector(model);
    compute_diffraction_intensity(model);
    unique();
    printf("Ending...\n");
}

XRD::~XRD()
{
    if(0!=numk){
        deallocate_2d(kvector, numk);
        deallocate(Ktheta);
        deallocate(Kintensity);
        deallocate(multiplicity);
    }
}

void XRD::count_diffraction_vector(MODEL *model)
{
    int c_kmin[3], c_kmax[3];
    vector_copy(c_kmin, kmin); vector_copy(c_kmax, kmax);
    for(int ih=c_kmin[0];ih<=c_kmax[0];ih++){
        for(int ik=c_kmin[1];ik<=c_kmax[1];ik++){
            for(int il=c_kmin[2];il<=c_kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double k[3]={double(ih), double(ik), double(il)};
                double Kmag=model->get_reciprocal_vector_length(k);
                if(2.0>=Kmag*model->lambda){
                    if(ih<kmin[0]) kmin[0]=ih;
                    if(ik<kmin[1]) kmin[1]=ik;
                    if(il<kmin[2]) kmin[2]=il;
                    if(ih>kmax[0]) kmax[0]=ih;
                    if(ik>kmax[1]) kmax[1]=ik;
                    if(il>kmax[2]) kmax[2]=il;
                    numk++;
                }
            }
        }
    }
    printf("Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
}

void XRD::compute_diffraction_intensity(MODEL *model)
{
    callocate_2d(&kvector, numk, 3, 0);
    callocate(&Ktheta, numk, 0.0);
    callocate(&Kintensity, numk, 0.0);
    printf("Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double k[3]={double(ih), double(ik), double(il)};
                double Kmag=model->get_reciprocal_vector_length(k);
                if(2>=Kmag*model->lambda){
                    double theta=asin(0.5*model->lambda*Kmag);
                    model->reciprocal_to_cartesian(k, k);
                    complex<double> factor=model->get_atomic_structure_factor(theta, k, is_lorentz_flag);
                    double intensity=model->get_diffraction_intensity(factor);
                    if(intensity>ZERO_LIMIT){
                        if(intensity<intensity_min) intensity_min=intensity;
                        if(intensity>intensity_max) intensity_max=intensity;
                        kvector[countk][0]=ih; kvector[countk][1]=ik; kvector[countk][2]=il;
                        Ktheta[countk]=theta; Kintensity[countk]=intensity;
                        countk++;
                    }
                }
            }
        }
    }
    numk=countk;
    printf("Ending computation of diffraction intensity\n");
    finish=clock();
    printf("Computation time [s]: %.8f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("Number of diffraction intensity: %d\n", numk);
    printf("Range of diffraction intensity: %.8f %.8f,", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf(" including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void XRD::count_diffraction_vector(MODEL *model, double spacingK[3], double minTheta, double maxTheta)
{
    printf("Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("Range of diffraction angle in reciprocal space (degree): %.8f %.8f\n", 2.0*minTheta/DEG_TO_RAD, 2.0*maxTheta/DEG_TO_RAD);
    int c_kmin[3]={kmin[0], kmin[1], kmin[2]}, c_kmax[3]={kmax[0], kmax[1], kmax[2]};
    for(int ih=c_kmin[0];ih<=c_kmax[0];ih++){
        for(int ik=c_kmin[1];ik<=c_kmax[1];ik++){
            for(int il=c_kmin[2];il<=c_kmax[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(2>=Kmag*model->lambda){
                    double theta=asin(0.5*model->lambda*Kmag);
                    if((theta<=maxTheta)&&(theta>=minTheta)){
                        if(ih<kmin[0]) kmin[0]=ih;
                        if(ik<kmin[1]) kmin[1]=ik;
                        if(il<kmin[2]) kmin[2]=il;
                        if(ih>kmax[0]) kmax[0]=ih;
                        if(ik>kmax[1]) kmax[1]=ik;
                        if(il>kmax[2]) kmax[2]=il;
                        numk++;
                    }
                }
            }
        }
    }
    printf("Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("Number of diffraction vector in reciprocal space: %d\n", numk);
}

void XRD::compute_diffraction_intensity(MODEL *model, double spacingK[3], double minTheta, double maxTheta)
{
    callocate_2d(&kvector, numk, 3, 0);
    callocate(&Ktheta, numk, 0.0);
    callocate(&Kintensity, numk, 0.0);
    printf("Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(2>=Kmag*model->lambda){
                    double theta=asin(0.5*model->lambda*Kmag);
                    if((theta<=maxTheta)&&(theta>=minTheta)){
                        complex<double> factor=model->get_atomic_structure_factor(theta, K, is_lorentz_flag);
                        double intensity=model->get_diffraction_intensity(factor);
                        if(intensity<intensity_min) intensity_min=intensity;
                        if(intensity>intensity_max) intensity_max=intensity;
                        kvector[countk][0]=ih; kvector[countk][1]=ik; kvector[countk][2]=il;
                        Ktheta[countk]=theta; Kintensity[countk]=intensity;
                        countk++;
                        if(0==countk%1000){
                            printf("Completed diffraction intensity %d of %d\n", countk, numk);
                        }
                    }
                }
            }
        }
    }
    printf("Ending computation of diffraction intensity\n");
    finish=clock();
    printf("Computation time [s]: %.8f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("Number of diffraction intensity: %d\n", countk);
    printf("Range of diffraction intensity: %.8f %.8f,", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf(" including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void XRD::quick_sort(int low, int high)
{
    int i=low, j=high;
    double theta=Ktheta[i];
    double intensity=Kintensity[i];
    int k[3]; vector_copy(k, kvector[i]);
    while(i<j){
        while(i<j&&Ktheta[j]>=theta){
            j--;
        }
        Ktheta[i]=Ktheta[j];
        Kintensity[i]=Kintensity[j];
        vector_copy(kvector[i], kvector[j]);
        while(i<j&&Ktheta[i]<=theta){
            i++;
        }
        Ktheta[j]=Ktheta[i];
        Kintensity[j]=Kintensity[i];
        vector_copy(kvector[j], kvector[i]);
    }
    Ktheta[i]=theta;
    Kintensity[i]=intensity;
    vector_copy(kvector[i], k);
    if(i-1>low){
        quick_sort(low, i-1);
    }
    if(i+1<high){
        quick_sort(i+1, high);
    }
}

void XRD::unique()
{
    quick_sort(0, numk-1);
    callocate(&multiplicity, numk, 0);
    int slow=0, fast=0;
    while(fast<numk){
        if(fabs(Ktheta[fast]-Ktheta[slow])>ZERO_LIMIT) {
            slow++;
            Ktheta[slow]=Ktheta[fast];
            Kintensity[slow]=Kintensity[fast];
            vector_copy(kvector[slow], kvector[fast]);
        }
        multiplicity[slow]++;
        fast++;
    }
    numk=slow+1;
    for(int i=0;i<numk;i++){
        Kintensity[i]*=(double)multiplicity[i];
        if(Kintensity[i]<intensity_min) intensity_min=Kintensity[i];
        if(Kintensity[i]>intensity_max) intensity_max=Kintensity[i];
    }
    printf("Number of unique diffraction intensity: %d", numk);
    printf("Range of diffraction intensity: %.8f %.8f,", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf(" including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void XRD::xrd(const char *xrd_path)
{
    FILE *fp=nullptr;
    fp=fopen(xrd_path,"w");
    fprintf(fp,"# h\tk\tl\t2theta\tintensity\tintensity_norm\tmultiplicity\n");
    double constn=100.0/intensity_max;
    for(int i=0;i<numk;i++){
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%d\n", 
                    kvector[i][0], kvector[i][1], kvector[i][2], 
                    2*Ktheta[i]/DEG_TO_RAD, Kintensity[i], constn*Kintensity[i], multiplicity[i]);
        fflush(fp);
    }
    fclose(fp);
    printf("Information for x-ray pattern stored in %s.\n", xrd_path);
}