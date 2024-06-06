#include "SND.h"

SND::SND(NMODEL *model, double min2Theta, double max2Theta, double spacing[3], bool is_spacing_auto, bool is_lorentz)
{
    printf("[INFO] Starting computation of selected-area neutron diffraction...\n");
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    minTheta=min2Theta/2.0*DEG_TO_RAD; maxTheta=max2Theta/2.0*DEG_TO_RAD;
    double Kmagnitude_max=2.0/model->lambda*sin(maxTheta);
    for(int i=0;i<3;i++){
        int NspacingK=ceil(Kmagnitude_max/spacingK[i]);
        kmin[i]=-NspacingK; kmax[i]=NspacingK;
    }
    is_lorentz_flag=is_lorentz;
    compute_diffraction_intensity(model);
    quick_unique();
    printf("[INFO] Ending computation of selected-area neutron diffraction\n");
}

SND::~SND()
{
    SND_KNODE *cur=khead;
    while(cur!=nullptr){
        SND_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void SND::add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity)
{
    if(ktail==nullptr){
        khead=ktail=new SND_KNODE;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;
        numk++;
    }else{
        ktail->next=new SND_KNODE;
        ktail=ktail->next;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;  
        numk++;
    }
}

void SND::compute_diffraction_intensity(NMODEL *model)
{
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Range of diffraction angle in reciprocal space (degree): %.8f %.8f\n", 2.0*minTheta/DEG_TO_RAD, 2.0*maxTheta/DEG_TO_RAD);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(2.0>=Kmag*model->lambda){
                    model->reciprocal_to_cartesian(K, K);
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, K, is_lorentz_flag);
                    if((intensity>ZERO_LIMIT)&&(theta<=maxTheta)&&(theta>=minTheta)){
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
    printf("[INFO] Ending computation of diffraction intensity\n");
    finish=clock();
    printf("[INFO] Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("[INFO] Number of diffraction intensity: %d\n", countk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f ", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf("including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void SND::copy_knode_data(SND_KNODE *knode1, SND_KNODE *knode2)
{
    vector_copy(knode1->k, knode2->k);
    knode1->theta=knode2->theta;
    knode1->intensity=knode2->intensity;
    knode1->multiplicity=knode2->multiplicity;
}

void SND::swap_knode_data(SND_KNODE *knode1, SND_KNODE *knode2)
{
    SND_KNODE *ktemp=new SND_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void SND::quick_sort(SND_KNODE *kstart, SND_KNODE *kend)
{
    if(kstart==nullptr||kend==nullptr||kstart==kend) return;
    SND_KNODE *knode1=kstart;
    SND_KNODE *knode2=kstart->next;
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

void SND::quick_unique()
{
    quick_sort(khead, ktail);
    int countk=0;
    SND_KNODE *kslow=khead, *kfast=khead;
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

    SND_KNODE *ktemp=khead;
    intensity_min=1.0e8; intensity_max=0.0;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        ktemp->intensity*=(double)ktemp->multiplicity;
        if(ktemp->intensity<intensity_min) intensity_min=ktemp->intensity;
        if(ktemp->intensity>intensity_max) intensity_max=ktemp->intensity;
        ktemp=ktemp->next;
    }
    printf("[INFO] Number of unique diffraction intensity: %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f ", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf("including lorentz factor\n");
    }else{
        printf("\n");
    }
}

void SND::snd(const char *snd_path, int nbin)
{
    FILE *fp=nullptr;
    fp=fopen(snd_path,"w");
    if(0==nbin){
        double constn=100.0/intensity_max;
        SND_KNODE *ktemp=khead;
        fprintf(fp,"# 2theta\tintensity\tintensity_norm (%d points)\n", numk);
        for(int i=0;i<numk&&ktemp!=nullptr;i++){
            fprintf(fp, "%.8f\t%.8f\t%.8f\n",
                        2*ktemp->theta/DEG_TO_RAD, ktemp->intensity, constn*ktemp->intensity);
            fflush(fp);
            ktemp=ktemp->next;
        }
    }else{
        double *kintensity=nullptr;
        callocate(&kintensity, nbin, 0.0);
        SND_KNODE *ktemp=khead;
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
        double constn=100.0/imax;
        fprintf(fp,"# 2theta\tintensity\tintensity_norm (%d bins)\n", nbin);
        for(int i=0;i<nbin;i++){
            fprintf(fp, "%.8f\t%.8f\t%.8f\n",
                        2*(minTheta+tbin*(i+0.5))/DEG_TO_RAD, kintensity[i], constn*kintensity[i]);
            fflush(fp);
        }
        deallocate(kintensity);
    }
    fclose(fp);
    printf("[INFO] Information for selected-area neutron pattern stored in %s\n", snd_path);
}