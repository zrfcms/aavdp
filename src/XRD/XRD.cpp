#include "XRD.h"

XRD::XRD(MODEL *model, double spacing[3], double min2Theta, double max2Theta, bool is_lorentz, bool is_spacing_auto)
{
    printf("Starting...\n");
    double spacingK[3];
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    minTheta=min2Theta/2.0*DEG_TO_RAD, maxTheta=max2Theta/2.0*DEG_TO_RAD;
    double Kmagnitude_max=2.0*sin(maxTheta)/model->lambda;
    for(int i=0;i<3;i++){
        int NspacingK=ceil(Kmagnitude_max/spacingK[i]);
        kmin[i]=-NspacingK; kmax[i]=NspacingK;
    }
    is_lorentz_flag=is_lorentz;
    compute_diffraction_intensity(model, spacingK);
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
    compute_diffraction_intensity(model);
    quick_unique();
    printf("Ending...\n");
}

XRD::~XRD()
{
    KNODE *cur=khead;
    while(cur!=nullptr){
        KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void XRD::add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity)
{
    if(ktail==nullptr){
        khead=ktail=new KNODE;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;
        numk++;
    }else{
        ktail->next=new KNODE;
        ktail=ktail->next;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;  
        numk++;
    }
}

void XRD::compute_diffraction_intensity(MODEL *model)
{
    printf("Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double k[3]={double(ih)*model->dimensionK[0], double(ik)*model->dimensionK[1], double(il)*model->dimensionK[2]};
                double Kmag=model->get_reciprocal_vector_length(k);
                if(2>=Kmag*model->lambda){
                    model->reciprocal_to_cartesian(k, k);
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, k, is_lorentz_flag);
                    if(intensity>ZERO_LIMIT){
                        if(ih<kmin[0]) kmin[0]=ih;
                        if(ik<kmin[1]) kmin[1]=ik;
                        if(il<kmin[2]) kmin[2]=il;
                        if(ih>kmax[0]) kmax[0]=ih;
                        if(ik>kmax[1]) kmax[1]=ik;
                        if(il>kmax[2]) kmax[2]=il;
                        if(intensity<intensity_min) intensity_min=intensity;
                        if(intensity>intensity_max) intensity_max=intensity;
                        add_k_node(ih, ik, il, theta, intensity, 0);
                        countk++;
                    }
                }
            }
        }
    }
    numk=countk;
    printf("Ending computation of diffraction intensity\n");
    finish=clock();
    printf("Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("Number of diffraction intensity: %d\n", numk);
    printf("Range of diffraction intensity: %.8f %.8f,", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf(" including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void XRD::compute_diffraction_intensity(MODEL *model, double spacingK[3])
{
    printf("Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("Range of diffraction angle in reciprocal space (degree): %.8f %.8f\n", 2.0*minTheta/DEG_TO_RAD, 2.0*maxTheta/DEG_TO_RAD);
    printf("Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
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
                        double intensity=model->get_diffraction_intensity(theta, K, is_lorentz_flag);
                        if(intensity<intensity_min) intensity_min=intensity;
                        if(intensity>intensity_max) intensity_max=intensity;
                        add_k_node(ih, ik, il, theta, intensity, 0);
                        countk++;
                    }
                }
            }
        }
    }
    numk=countk;
    printf("Ending computation of diffraction intensity\n");
    finish=clock();
    printf("Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("Number of diffraction intensity: %d\n", countk);
    printf("Range of diffraction intensity: %.8f %.8f,", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf(" including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void XRD::copy_knode_data(KNODE *knode1, KNODE *knode2)
{
    vector_copy(knode1->k, knode2->k);
    knode1->theta=knode2->theta;
    knode1->intensity=knode2->intensity;
    knode1->multiplicity=knode2->multiplicity;
}

void XRD::swap_knode_data(KNODE *knode1, KNODE *knode2)
{
    KNODE *ktemp=new KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void XRD::quick_sort(KNODE *kstart, KNODE *kend)
{
    if(kstart==nullptr||kend==nullptr||kstart==kend) return;
    KNODE *knode1=kstart;
    KNODE *knode2=kstart->next;
    double theta=kstart->theta;
    while(knode2!=kend->next&&knode2!=nullptr){
        if(knode2->theta<theta){
            knode1=knode1->next;
            if(knode1!=knode2){
                swap_knode_data(knode1, knode2);
            }
        }
        knode2=knode2->next;
    }
    swap_knode_data(knode1, kstart);
    quick_sort(kstart, knode1);
    quick_sort(knode1->next, kend);
}

void XRD::quick_unique()
{
    quick_sort(khead, ktail);
    int countk=0;
    KNODE *kslow=khead, *kfast=khead;
    while(kfast!=nullptr){
        if(fabs(kfast->theta-kslow->theta)>ZERO_LIMIT) {
            kslow->next=kfast;
            kslow=kslow->next;
            countk++;
        }
        kslow->multiplicity++;
        kfast=kfast->next;
    }
    numk=countk+1;

    KNODE *ktemp=khead;
    intensity_min=1.0e8; intensity_max=0.0;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        ktemp->intensity*=(double)ktemp->multiplicity;
        if(ktemp->intensity<intensity_min) intensity_min=ktemp->intensity;
        if(ktemp->intensity>intensity_max) intensity_max=ktemp->intensity;
        ktemp=ktemp->next;
    }
    printf("Number of unique diffraction intensity: %d\n", numk);
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
    fprintf(fp,"# 2theta\tintensity\tintensity_norm\n");
    double constn=100.0/intensity_max;
    KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        fprintf(fp, "%.8f\t%.8f\t%.8f\n",
                    2*ktemp->theta/DEG_TO_RAD, ktemp->intensity, constn*ktemp->intensity);
        fflush(fp);
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("Information for x-ray pattern stored in %s.\n", xrd_path);
}

void XRD::xrd(const char *xrd_path, int nbin)
{
    if(0==nbin){
        xrd(xrd_path); return;
    }
    double *kintensity=nullptr;
    callocate(&kintensity, nbin, 0.0);
    KNODE *ktemp=khead;
    double tbin=(maxTheta-minTheta)/double(nbin);
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        int j=int((ktemp->theta-minTheta)/tbin);
        kintensity[j]+=ktemp->intensity;
        ktemp=ktemp->next;
    }
    double imax=kintensity[0];
    for(int i=1;i<nbin;i++){
        if(imax<kintensity[i]) imax=kintensity[i];
    }

    FILE *fp=nullptr;
    fp=fopen(xrd_path,"w");
    fprintf(fp,"# 2theta\tintensity\tintensity_norm\n");
    double constn=100.0/imax;
    for(int i=0;i<nbin;i++){
        fprintf(fp, "%.8f\t%.8f\t%.8f\n",
                    2*(minTheta+tbin*(i+0.5))/DEG_TO_RAD, kintensity[i], constn*kintensity[i]);
        fflush(fp);
    }
    fclose(fp);
    deallocate(kintensity);
    printf("Information for x-ray pattern stored in %s.\n", xrd_path);
}

// void XRD::quick_sort(int low, int high)
// {
//     int i=low, j=high;
//     double theta=Ktheta[i];
//     double intensity=Kintensity[i];
//     int k[3]; vector_copy(k, kvector[i]);
//     while(i<j){
//         while(i<j&&Ktheta[j]>=theta){
//             j--;
//         }
//         Ktheta[i]=Ktheta[j];
//         Kintensity[i]=Kintensity[j];
//         vector_copy(kvector[i], kvector[j]);
//         while(i<j&&Ktheta[i]<=theta){
//             i++;
//         }
//         Ktheta[j]=Ktheta[i];
//         Kintensity[j]=Kintensity[i];
//         vector_copy(kvector[j], kvector[i]);
//     }
//     Ktheta[i]=theta;
//     Kintensity[i]=intensity;
//     vector_copy(kvector[i], k);
//     if(i-1>low){
//         quick_sort(low, i-1);
//     }
//     if(i+1<high){
//         quick_sort(i+1, high);
//     }
// }

// void XRD::unique()
// {
//     quick_sort(0, numk-1);
//     callocate(&multiplicity, numk, 0);
//     int slow=0, fast=0;
//     while(fast<numk){
//         if(fabs(Ktheta[fast]-Ktheta[slow])>ZERO_LIMIT) {
//             slow++;
//             Ktheta[slow]=Ktheta[fast];
//             Kintensity[slow]=Kintensity[fast];
//             vector_copy(kvector[slow], kvector[fast]);
//         }
//         multiplicity[slow]++;
//         fast++;
//     }
//     numk=slow+1;
//     for(int i=0;i<numk;i++){
//         Kintensity[i]*=(double)multiplicity[i];
//         if(Kintensity[i]<intensity_min) intensity_min=Kintensity[i];
//         if(Kintensity[i]>intensity_max) intensity_max=Kintensity[i];
//     }
//     printf("Number of unique diffraction intensity: %d", numk);
//     printf("Range of diffraction intensity: %.8f %.8f,", intensity_min, intensity_max);
//     if(is_lorentz_flag){
//         printf(" including lorentz factor\n");
//     }else{
//         printf("\n");
//     }
// }