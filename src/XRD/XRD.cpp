#include "XRD.h"

XRD::XRD(XMODEL *model, double min2Theta, double max2Theta, int lp_type, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of x-ray diffraction...\n");
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    minTheta=min2Theta/2.0*DEG_TO_RAD; maxTheta=max2Theta/2.0*DEG_TO_RAD;
    double Kmagnitude_max=2.0/model->lambda*sin(maxTheta);
    for(int i=0;i<3;i++){
        int NspacingK=ceil(Kmagnitude_max/spacingK[i]);
        kmin[i]=-NspacingK; kmax[i]=NspacingK;
    }
    compute_diffraction_intensity(model, lp_type);
    quick_unique();
    printf("[INFO] Ending computation of x-ray diffraction\n");
}

XRD::~XRD()
{
    XRD_KNODE *cur=khead;
    while(cur!=nullptr){
        XRD_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void XRD::add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity)
{
    if(ktail==nullptr){
        khead=ktail=new XRD_KNODE;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;
        numk++;
    }else{
        ktail->next=new XRD_KNODE;
        ktail=ktail->next;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;  
        numk++;
    }
}

void XRD::compute_diffraction_intensity(XMODEL *model, int lp_type)
{
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int num=(kmax[0]-kmin[0]+1)*(kmax[1]-kmin[1]+1)*(kmax[2]-kmin[2]+1);
    int count=0;
    int kkmin[3], kkmax[3];
    vector_copy(kkmin, kmin);
    vector_copy(kkmax, kmax);
    for(int ih=kkmin[0];ih<=kkmax[0];ih++){
        for(int ik=kkmin[1];ik<=kkmax[1];ik++){
            for(int il=kkmin[2];il<=kkmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(2.0>=Kmag*model->lambda){
                    model->reciprocal_to_cartesian(K, K);
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, K, lp_type);
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
                    }
                }
                count++;
                if(0==count%100){
                    printf("[INFO] Completed diffraction intensity %d of %d\n", count, num);
                }
            }
        }
    }
    printf("[INFO] Ending computation of diffraction intensity\n");
    finish=clock();
    printf("[INFO] Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("[INFO] Number of diffraction intensity: %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void XRD::copy_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2)
{
    vector_copy(knode1->k, knode2->k);
    knode1->theta=knode2->theta;
    knode1->intensity=knode2->intensity;
    knode1->multiplicity=knode2->multiplicity;
}

void XRD::swap_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2)
{
    XRD_KNODE *ktemp=new XRD_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void XRD::quick_sort(XRD_KNODE *kstart, XRD_KNODE *kend)
{
    if(kstart==nullptr||kend==nullptr||kstart==kend) return;
    XRD_KNODE *knode1=kstart;
    XRD_KNODE *knode2=kstart->next;
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
    XRD_KNODE *kslow=khead, *kfast=khead;
    while(kfast!=nullptr){
        if(fabs(kfast->theta-kslow->theta)>ZERO_LIMIT){
            kslow->next=kfast;
            kslow=kslow->next;
            countk++;
        }
        kslow->multiplicity++;
        kfast=kfast->next;
    }
    numk=countk+1;

    XRD_KNODE *ktemp=khead;
    intensity_min=1.0e8; intensity_max=0.0;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        ktemp->intensity*=(double)ktemp->multiplicity;
        if(ktemp->intensity<intensity_min) intensity_min=ktemp->intensity;
        if(ktemp->intensity>intensity_max) intensity_max=ktemp->intensity;
        ktemp=ktemp->next;
    }
    printf("[INFO] Number of unique diffraction intensity: %d\n", numk);
    printf("[INFO] Range of unique diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void XRD::xrd(char *xrd_path, char *png_path)
{
    FILE *fp=nullptr;
    fp=fopen(xrd_path,"w");
    double constn=100.0/intensity_max;
    fprintf(fp,"# H\tK\tL\tmultiplicity\t2theta\tintensity\tintensity_norm (%d points)\n", numk);
    double *ktheta=nullptr, *kintensity=nullptr;
    callocate(&ktheta, numk, 0.0);
    callocate(&kintensity, numk, 0.0);
    XRD_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        ktheta[i]=ktemp->theta*2.0*RAD_TO_DEG; kintensity[i]=constn*ktemp->intensity;
        fprintf(fp, "%d\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\n", ktemp->k[0], ktemp->k[1], ktemp->k[2], 
                ktemp->multiplicity, ktheta[i], ktemp->intensity, kintensity[i]);
        fflush(fp);
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for x-ray pattern stored in %s\n", xrd_path);

    double min2Theta=2*minTheta/DEG_TO_RAD, max2Theta=2*maxTheta/DEG_TO_RAD;
    int n_major_xtick=19, n_major_ytick=11;
    double *major_xticks, *major_yticks;
    mallocate(&major_xticks, n_major_xtick);
    mallocate(&major_yticks, n_major_ytick);
    for(int i=0;i<=18;i++){
        major_xticks[i]=double(i)*10.0;
    }
    for(int i=0;i<=10;i++){
        major_yticks[i]=double(i)*10.0;
    }
    double height=4.0, width=6.0;
    GRAPH graph(width, height, 300);
    graph.set_xlim(min2Theta, max2Theta);
    graph.set_ylim(0.0, 100.0);
    graph.set_xticks(major_xticks, n_major_xtick);
    graph.set_yticks(major_yticks, n_major_ytick);
    graph.set_tick_in(false);
    graph.set_top_visible(false);
    graph.set_right_visible(false);
    graph.hist(ktheta, kintensity, numk);
    graph.draw(png_path);
    deallocate(ktheta);
    deallocate(kintensity);
    printf("[INFO] Image for x-ray pattern stored in %s\n", png_path);
}

void XRD::xrd(char *xrd_path, char *png_path, double peak_parameter, double scherrer_lambda, double scherrer_size, double dt)
{
    double tbin=dt/2.0*DEG_TO_RAD;
    int    nbin=round((maxTheta-minTheta)/tbin);
    double *ktheta=nullptr, *kintensity=nullptr;
    callocate(&ktheta, nbin, 0.0);
    callocate(&kintensity, nbin, 0.0);
    for(int i=0;i<nbin;i++){
        ktheta[i]=(minTheta+tbin*(i+0.5));
    }

    double constw=SCHERRER_CONST*scherrer_lambda/scherrer_size;
    double *kweights=nullptr, *kintensity_c=nullptr; 
    callocate(&kweights, nbin, 0.0);
    callocate(&kintensity_c, nbin, 0.0);
    XRD_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        int j=floor((ktemp->theta-minTheta)/tbin);
        double FWHM=constw/cos(ktheta[j]);
        kintensity_c[j]=ktemp->intensity;
        pseudo_Voigt(kweights, ktheta, nbin, peak_parameter, ktheta[j], FWHM);
        convolve(kintensity_c, kintensity_c, kweights, nbin, nbin, j);
        for(int k=0;k<nbin;k++){
            kintensity[k]+=kintensity_c[k];
            kintensity_c[k]=0.0; kweights[k]=0.0;
        }
        ktemp=ktemp->next;
    }
    deallocate(kweights);
    deallocate(kintensity_c);

    double imax=0.0, imin=1.0e8;
    for(int i=0;i<nbin;i++){
        if(imax<kintensity[i]) imax=kintensity[i];
        if(imin>kintensity[i]) imin=kintensity[i];
    }
    printf("[INFO] Number of profiled diffraction intensity: %d\n", nbin);
    printf("[INFO] Range of profiled diffraction intensity: %.8f %.8f\n", imin, imax);
    double constn=100.0/imax;
    FILE *fp=nullptr;
    fp=fopen(xrd_path,"w");
    fprintf(fp,"# 2theta\tintensity\tintensity_norm (%d bins)\n", nbin);
    for(int i=0;i<nbin;i++){
        ktheta[i]*=(2.0*RAD_TO_DEG);
        fprintf(fp, "%.8f\t%.8f\t", ktheta[i], kintensity[i]);
        kintensity[i]*=constn;
        fprintf(fp, "%.8f\n", kintensity[i]);
        fflush(fp);
    }
    fclose(fp);
    printf("[INFO] Information for x-ray pattern stored in %s\n", xrd_path);

    double min2Theta=2.0*minTheta/DEG_TO_RAD, max2Theta=2.0*maxTheta/DEG_TO_RAD;
    int n_major_xtick=19, n_major_ytick=11;
    double *major_xticks, *major_yticks;
    mallocate(&major_xticks, n_major_xtick);
    mallocate(&major_yticks, n_major_ytick);
    for(int i=0;i<=18;i++){
        major_xticks[i]=double(i)*10.0;
    }
    for(int i=0;i<=10;i++){
        major_yticks[i]=double(i)*10.0;
    }
    double height=4.0, width=6.0;
    GRAPH graph(width, height, 300);
    graph.set_xlim(min2Theta, max2Theta);
    graph.set_ylim(0.0, 100.0);
    graph.set_xticks(major_xticks, n_major_xtick);
    graph.set_yticks(major_yticks, n_major_ytick);
    graph.set_tick_in(false);
    graph.set_top_visible(false);
    graph.set_right_visible(false);
    graph.line(ktheta, kintensity, nbin);
    graph.draw(png_path);
    deallocate(ktheta);
    deallocate(kintensity);
    printf("[INFO] Image for x-ray pattern stored in %s\n", png_path);
}