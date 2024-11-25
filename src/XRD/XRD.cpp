#include "XRD.h"

void copy_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2)
{
    vector_copy(knode1->hkl, knode2->hkl);
    knode1->theta=knode2->theta;
    knode1->intensity=knode2->intensity;
    knode1->multiplicity=knode2->multiplicity;
}

void swap_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2)
{
    XRD_KNODE *ktemp=new XRD_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void quick_sort(XRD_KNODE *khead, XRD_KNODE *kend)
{
    if(khead==nullptr||kend==nullptr||khead==kend) return;
    XRD_KNODE *knode1=khead;
    XRD_KNODE *knode2=khead->next;
    double theta=khead->theta;
    while(knode2!=kend->next&&knode2!=nullptr){
        if(knode2->theta<theta){
            knode1=knode1->next;
            if(knode1!=knode2){
                swap_knode_data(knode1, knode2);
            }
        }
        knode2=knode2->next;
    }
    swap_knode_data(knode1, khead);
    quick_sort(khead, knode1);
    quick_sort(knode1->next, kend);
}

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

XRD::XRD(MODEL *model, double min2Theta, double max2Theta, double threshold, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of %s diffraction...\n", model->radiation);
    minTheta=min2Theta*DEG_TO_RAD_HALF; 
    maxTheta=max2Theta*DEG_TO_RAD_HALF;
    compute_diffraction_intensity(model, spacing, is_spacing_auto);
    filter_diffraction_intensity(threshold);
    printf("[INFO] Number of diffraction intensity: %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    unique_diffraction_intensity();
    printf("[INFO] Number of unique diffraction intensity: %d\n", numk);
    printf("[INFO] Range of unique diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    printf("[INFO] Ending computation of %s diffraction\n", model->radiation);
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

void XRD::add_k_node(int hkl[3], double theta, double intensity, int multiplicity)
{
    if(ktail==nullptr){
        khead=ktail=new XRD_KNODE;
        vector_copy(ktail->hkl, hkl);
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;
        numk++;
    }else{
        ktail->next=new XRD_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        ktail->theta=theta; ktail->intensity=intensity;
        ktail->multiplicity=multiplicity;  
        numk++;
    }
    if(intensity<intensity_min) intensity_min=intensity;
    if(intensity>intensity_max) intensity_max=intensity;
}

void XRD::compute_diffraction_intensity(MODEL *model, double spacing[3], bool is_spacing_auto)
{
    double Kmagnitude_min=2.0/model->lambda*sin(minTheta);
    double Kmagnitude_max=2.0/model->lambda*sin(maxTheta);
    double spacingK[3];
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    printf("[INFO] Range along a*, b*, or c* in reciprocal space (Angstrom-1): %.8f %.8f\n", Kmagnitude_min, Kmagnitude_max);
    printf("[INFO] Spacings along a*, b*, and c* in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    int kmin[3], kmax[3];
    for(int i=0;i<3;i++){
        kmin[i]=floor(Kmagnitude_min/spacingK[i]); 
        kmax[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    int num=(2*kmax[0]+1)*(2*kmax[1]+1)*(2*kmax[2]+1);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int count=0;
    for(int ih=-kmax[0];ih<=kmax[0];++ih){
        for(int ik=-kmax[1];ik<=kmax[1];++ik){
            for(int il=-kmax[2];il<=kmax[2];++il){
                if(0==ih&&0==ik&&0==il) continue;
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(Kmag<=Kmagnitude_max&&Kmag>=Kmagnitude_min){
                    model->reciprocal_to_cartesian(K, K);
                    int    hkl[3]={ih, ik, il};
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, K);
                    if(intensity>XRD_INTENSITY_LIMIT){
                        add_k_node(hkl, theta, intensity, 1);
                    }
                }
                count++;
                if(0==count%1000){
                    printf("[INFO] Completed diffraction intensity %d of %d\n", count, num);
                }
            }
        }
    }
    finish=clock();
    printf("[INFO] Ending computation of diffraction intensity\n");
    printf("[INFO] Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);
}

void XRD::filter_diffraction_intensity(double threshold)
{
    intensity_min=1.0e8;
    double intensity_threshold=intensity_max*threshold;
    while(khead!=nullptr&&khead->intensity<=intensity_threshold){
        XRD_KNODE *ktemp=khead;
        khead=khead->next;
        delete ktemp;
        numk--;
    }
    XRD_KNODE *cur=khead;
    while(cur!=nullptr&&cur->next!=nullptr){
        if(cur->next->intensity<=intensity_threshold){
            XRD_KNODE *ktemp=cur->next;
            cur->next=cur->next->next;
            delete ktemp;
            numk--;
        }else{
            if(intensity_min>cur->intensity) intensity_min=cur->intensity;
            cur=cur->next;
        }
    }
}

void XRD::unique_diffraction_intensity()
{
    quick_sort(khead, ktail);
    int countk=1;
    XRD_KNODE *kslow=khead, *kfast=khead->next;
    while(kfast!=nullptr){
        if(fabs(kfast->theta-kslow->theta)>XRD_THETA_LIMIT){
            kslow->next=kfast;
            kslow=kslow->next;
            kfast=kfast->next;
            countk++;
        }else{
            kslow->multiplicity++;
            kslow->next=nullptr;
            XRD_KNODE *ktemp=kfast;
            kfast=kfast->next;
            free(ktemp);
        }
    }
    numk=countk;

    XRD_KNODE *ktemp=khead;
    intensity_min=1.0e8; intensity_max=0.0;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        ktemp->intensity*=(double)ktemp->multiplicity;
        if(ktemp->intensity<intensity_min) intensity_min=ktemp->intensity;
        if(ktemp->intensity>intensity_max) intensity_max=ktemp->intensity;
        ktemp=ktemp->next;
    }
}

void XRD::img(char *png_path, double *x, double *y, int num, double xmin, double xmax, char mode)
{
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
    graph.set_xlim(xmin, xmax);
    graph.set_ylim(0.0, 100.0);
    graph.set_xticks(major_xticks, n_major_xtick);
    graph.set_yticks(major_yticks, n_major_ytick);
    graph.set_tick_in(false);
    graph.set_top_visible(false);
    graph.set_right_visible(false);
    switch(mode){
        case 'h':
            graph.hist(x, y, num);
            break;
        case 'l':
            graph.line(x, y, num);
            break;
        default:
            graph.hist(x, y, num);
            printf("[WARNING] Unrecognized image mode '%c'", mode);
    }
    graph.draw(png_path);
}

void XRD::xrd(char *xrd_path)
{
    FILE *fp=fopen(xrd_path,"w");
    fprintf(fp,"# h\tk\tl\tmultiplicity\t2theta\tintensity\tintensity_norm (%d points)\n", numk);
    XRD_KNODE *ktemp=khead;
    double constn=100.0/intensity_max;
    double *theta=nullptr, *intensity=nullptr;
    callocate(&theta, numk, 0.0); callocate(&intensity, numk, 0.0);
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        theta[i]=ktemp->theta*RAD_TO_DEG_TWO;
        intensity[i]=ktemp->intensity*constn;
        fprintf(fp, "%d\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\n", ktemp->hkl[0], ktemp->hkl[1], ktemp->hkl[2], ktemp->multiplicity, theta[i], ktemp->intensity, intensity[i]);
        fflush(fp);
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", xrd_path);

    char png_path[PATH_CHAR_NUMBER];
    strcpy(png_path, xrd_path); strcat(png_path, ".png");
    img(png_path, theta, intensity, numk, minTheta*RAD_TO_DEG_TWO, maxTheta*RAD_TO_DEG_TWO, 'h');
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void XRD::xrd(char *xrd_path, double mixing_param, double scherrer_lambda, double scherrer_diameter, double bin2Theta)
{
    double tbin=bin2Theta*DEG_TO_RAD_HALF;
    int    nbin=round((maxTheta-minTheta)/tbin);
    double *theta=nullptr, *intensity=nullptr;
    callocate(&theta, nbin, 0.0);
    callocate(&intensity, nbin, 0.0);
    for(int i=0;i<nbin;i++){
        theta[i]=(minTheta+tbin*(i+0.5));
    }
    double constw=SCHERRER_CONST*scherrer_lambda/scherrer_diameter;
    double *intensity_c=nullptr; 
    callocate(&intensity_c, nbin, 0.0);
    XRD_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        int j=floor((ktemp->theta-minTheta)/tbin);
        double FWHM=constw/cos(theta[j]);
        pseudo_Voigt(intensity_c, theta, nbin, theta[j], ktemp->intensity, mixing_param, FWHM);
        for(int k=0;k<nbin;k++){
            intensity[k]+=intensity_c[k];
            intensity_c[k]=0.0;
        }
        ktemp=ktemp->next;
    }
    deallocate(intensity_c);
    double imax=0.0, imin=1.0e8;
    for(int i=0;i<nbin;i++){
        if(imax<intensity[i]) imax=intensity[i];
        if(imin>intensity[i]) imin=intensity[i];
    }
    printf("[INFO] Number of profiled diffraction intensity: %d\n", nbin);
    printf("[INFO] Range of profiled diffraction intensity: %.8f %.8f\n", imin, imax);

    FILE *fp=nullptr;
    fp=fopen(xrd_path,"w");
    fprintf(fp,"# 2theta\tintensity\tintensity_norm (%d bins)\n", nbin);
    double constn=100.0/imax;
    for(int i=0;i<nbin;i++){
        theta[i]*=RAD_TO_DEG_TWO;
        fprintf(fp, "%.8f\t%.8f\t", theta[i], intensity[i]);
        intensity[i]*=constn;
        fprintf(fp, "%.8f\n", intensity[i]);
        fflush(fp);
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", xrd_path);

    char png_path[PATH_CHAR_NUMBER];
    strcpy(png_path, xrd_path); strcat(png_path, ".png");
    img(png_path, theta, intensity, nbin, minTheta*RAD_TO_DEG_TWO, maxTheta*RAD_TO_DEG_TWO, 'l');
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void XRD::xrd(char *xrd_path, double mixing_param, double FWHM, double bin2Theta)
{
    double tbin=bin2Theta*DEG_TO_RAD_HALF;
    int    nbin=round((maxTheta-minTheta)/tbin);
    double *theta=nullptr, *intensity=nullptr;
    callocate(&theta, nbin, 0.0);
    callocate(&intensity, nbin, 0.0);
    for(int i=0;i<nbin;i++){
        theta[i]=(minTheta+tbin*(i+0.5));
    }
    double *intensity_c=nullptr; 
    callocate(&intensity_c, nbin, 0.0);
    XRD_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        int j=floor((ktemp->theta-minTheta)/tbin);
        pseudo_Voigt(intensity_c, theta, nbin, theta[j], ktemp->intensity, mixing_param, FWHM*DEG_TO_RAD);
        for(int k=0;k<nbin;k++){
            intensity[k]+=intensity_c[k];
            intensity_c[k]=0.0;
        }
        ktemp=ktemp->next;
    }
    deallocate(intensity_c);
    double imax=0.0, imin=1.0e8;
    for(int i=0;i<nbin;i++){
        if(imax<intensity[i]) imax=intensity[i];
        if(imin>intensity[i]) imin=intensity[i];
    }
    printf("[INFO] Number of profiled diffraction intensity: %d\n", nbin);
    printf("[INFO] Range of profiled diffraction intensity: %.8f %.8f\n", imin, imax);

    FILE *fp=nullptr;
    fp=fopen(xrd_path,"w");
    fprintf(fp,"# 2theta\tintensity\tintensity_norm (%d bins)\n", nbin);
    double constn=100.0/imax;
    for(int i=0;i<nbin;i++){
        theta[i]*=RAD_TO_DEG_TWO;
        fprintf(fp, "%.8f\t%.8f\t", theta[i], intensity[i]);
        intensity[i]*=constn;
        fprintf(fp, "%.8f\n", intensity[i]);
        fflush(fp);
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", xrd_path);

    char png_path[PATH_CHAR_NUMBER];
    strcpy(png_path, xrd_path); strcat(png_path, ".png");
    img(png_path, theta, intensity, nbin, minTheta*RAD_TO_DEG_TWO, maxTheta*RAD_TO_DEG_TWO, 'l');
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

