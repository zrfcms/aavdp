#include "KED.h"

void copy_knode_data(KED_KNODE *knode1, KED_KNODE *knode2)
{
    vector_copy(knode1->hkl, knode2->hkl);
    vector_copy(knode1->K, knode2->K);
    knode1->Kmagnitude=knode2->Kmagnitude;
    knode1->intensity=knode2->intensity;
}

void swap_knode_data(KED_KNODE *knode1, KED_KNODE *knode2)
{
    KED_KNODE *ktemp=new KED_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void quick_sort(KED_KNODE *kstart, KED_KNODE *kend)
{
    if(kstart==nullptr||kend==nullptr||kstart==kend) return;
    KED_KNODE *knode1=kstart;
    KED_KNODE *knode2=kstart->next;
    double Kmagnitude=kstart->Kmagnitude;
    while(knode2!=kend->next&&knode2!=nullptr){
        if(knode2->Kmagnitude<Kmagnitude){
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

void quick_sort(KED_KNODE *kstart)
{
    if(kstart==nullptr) return;
    KED_KNODE *knode1=kstart;
    KED_KNODE *knode2=kstart->next;
    double hkl[3]; vector_copy(hkl, kstart->hkl);
    while(knode2!=nullptr){
        if((0==int(knode2->hkl[0])+int(hkl[0]))&&(0==int(knode2->hkl[1])+int(hkl[1]))&&(0==int(knode2->hkl[2])+int(hkl[2]))){
            knode1=knode1->next;
            if(knode1!=knode2){
                swap_knode_data(knode1, knode2);
            }
            break;
        }else{
            knode2=knode2->next;
        }
        
    }
    quick_sort(knode1->next);
}

KED::KED(MODEL *model, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of kinematical electron diffraction...\n");
    printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", model->lambda);
    lambda=model->lambda;
    Kmagnitude_max=Kmag_max;
    compute_diffraction_intensity(model, spacing, is_spacing_auto);
    printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
    printf("[INFO] Number of diffraction intensity (including intensity at the transmission spot): %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    filter_diffraction_intensity(threshold);
    quick_sort(khead, ktail);
    printf("[INFO] Number of filtered diffraction intensity (including intensity at the transmission spot): %d\n", numk);
    printf("[INFO] Range of filtered diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    printf("[INFO] Ending computation of kinematical electron diffraction\n");
}

KED::KED(MODEL *model, int zone[3], double thickness, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of kinematical electron diffraction...\n");
    printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", model->lambda);
    lambda=model->lambda;
    Kmagnitude_max=Kmag_max;
    compute_diffraction_intensity(model, zone, thickness, spacing, is_spacing_auto);
    printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
    printf("[INFO] Number of diffraction intensity (including intensity at the transmission spot): %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    filter_diffraction_intensity(threshold);
    printf("[INFO] Number of filtered diffraction intensity (including intensity at the transmission spot): %d\n", numk);
    printf("[INFO] Range of filtered diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    find_first_and_second_knearests();
    if(knearest_2==nullptr){
        printf("[INFO] The first nearest diffraction vectors R1: [%.8f %.8f %.8f]\n", knearest_1->K[0], knearest_1->K[1], knearest_1->K[2]);
        printf("[WARN] Unable to find the second nearest diffraction vector\n");
    }else{
        printf("[INFO] The first and second nearest diffraction vectors R1, R2: [%.8f %.8f %.8f], [%.8f %.8f %.8f] (R2/R1 %.8f and angle %.8f)\n", 
        knearest_1->K[0], knearest_1->K[1], knearest_1->K[2], knearest_2->K[0], knearest_2->K[1], knearest_2->K[2], knearest_2->Kmagnitude/knearest_1->Kmagnitude, vector_angle(knearest_2->K, knearest_1->K)*RAD_TO_DEG);
    }
    rotate_by_first_knearest(zone);
    printf("[INFO] Ending computation of kinematical electron diffraction\n");
}

KED::~KED()
{
    KED_KNODE *cur=khead;
    while(cur!=nullptr){
        KED_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void KED::add_k_node(double hkl[3], double K[3], double Kmagnitude, double intensity)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new KED_KNODE;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kmagnitude=Kmagnitude;
        ktail->intensity=intensity;
        numk++;
    }else{
        ktail->next=new KED_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kmagnitude=Kmagnitude;
        ktail->intensity=intensity;
        if(intensity_max<intensity) intensity_max=intensity;
        if(intensity_min>intensity) intensity_min=intensity;
        numk++;
    }
}

void KED::compute_diffraction_intensity(MODEL *model, double spacing[3], bool is_spacing_auto)
{
    double spacingK[3];
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    int kmin[3], kmax[3];
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    printf("[INFO] Spacings along a*, b*, and c* in reciprocal space [Angstrom-1]: %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Number of spacings along a*, b*, and c* in reciprocal space: %d %d %d\n", kmax[0], kmax[1], kmax[2]);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    double hkl[3];
    double K[3]={0.0};
    double intensity=model->get_diffraction_intensity(0.0, K, true);
    add_k_node(hkl, K, 0.0, intensity);
    clock_t start, finish;
    start=clock();
    int num=(2*kmax[0]+1)*(2*kmax[1]+1)*(2*kmax[2]+1);
    int count=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                count++;
                if(0==ih&&0==ik&&0==il) continue;
                double hkl[3]={double(ih), double(ik), double(il)};
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmagnitude=model->get_reciprocal_vector_length(K);
                if(Kmagnitude<Kmagnitude_max){
                    model->reciprocal_to_cartesian(K, K);
                    double intensity=model->get_diffraction_intensity(Kmagnitude, K, false);
                    if(intensity>KED_INTENSITY_LIMIT){
                        add_k_node(hkl, K, Kmagnitude, intensity);
                    }
                }
                if(0==count%1000){
                    printf("[INFO] Completed diffraction intensity %d of %d\n", count, num);
                }
            }
        }
    }
    printf("[INFO] Ending computation of diffraction intensity\n");
    finish=clock();
    printf("[INFO] Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);
}

void KED::compute_diffraction_intensity(MODEL *model, int zone[3], double thickness, double spacing[3], bool is_spacing_auto)
{
    double spacingK[3];
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    int kmin[3], kmax[3];
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    double upper_bound=thickness/2.0;
    double lower_bound=-thickness/2.0;
    printf("[INFO] Spacings along a*, b*, and c* in reciprocal space [Angstrom-1]: %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Number of spacings along a*, b*, and c* in reciprocal space: %d %d %d\n", kmax[0], kmax[1], kmax[2]);
    printf("[INFO] Range along zone-[%.8f %.8f %.8f] in reciprocal space [Angstrom-1]: %.8f %.8f\n", n_zone[0], n_zone[1], n_zone[2], lower_bound, upper_bound);
    double hkl[3]={0.0};
    double K[3]={0.0};
    double intensity=model->get_diffraction_intensity(0.0, K, true);
    add_k_node(hkl, K, 0.0, intensity);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int num=(2*kmax[0]+1)*(2*kmax[1]+1)*(2*kmax[2]+1);
    int count=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double hkl[3]={double(ih), double(ik), double(il)};
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmagnitude=model->get_reciprocal_vector_length(K);
                if(Kmagnitude<Kmagnitude_max){
                    model->reciprocal_to_cartesian(K, K);
                    double proj=vector_dot(K, n_zone);
                    if((proj>lower_bound)&&(proj<upper_bound)){
                        double intensity=model->get_diffraction_intensity(Kmagnitude, K, false);
                        if(intensity>KED_INTENSITY_LIMIT){
                            add_k_node(hkl, K, Kmagnitude, intensity);
                        }
                    }
                }
                count++;
                if(0==count%1000){
                    printf("[INFO] Completed diffraction intensity %d of %d\n", count, num);
                }
            }
        }
    }
    printf("[INFO] Ending computation of diffraction intensity\n");
    finish=clock();
    printf("[INFO] Computation time [s]: %.2f.\n", double(finish-start)/CLOCKS_PER_SEC);
}

void KED::filter_diffraction_intensity(double threshold)
{
    double intensity_threshold=intensity_max*threshold;
    intensity_min=1.0e8;
    while(khead!=nullptr&&khead->intensity<=intensity_threshold){
        KED_KNODE *ktemp=khead;
        khead=khead->next;
        delete ktemp;
        numk--;
    }
    KED_KNODE *cur=khead;
    while(cur!=nullptr&&cur->next!=nullptr){
        if(cur->next->intensity<=intensity_threshold){
            KED_KNODE *ktemp=cur->next;
            cur->next=cur->next->next;
            delete ktemp;
            numk--;
        }else{
            if(intensity_min>cur->intensity) intensity_min=cur->intensity;
            cur=cur->next;
        }
    }
}

void KED::find_first_and_second_knearests()
{
    quick_sort(khead, ktail);
    knearest_1=khead->next;
    if(knearest_1==nullptr){
        printf("[ERROR] Unable to find the first nearest diffraction vector\n");
        exit(1);
    }
    knearest_2=knearest_1->next;
    while(knearest_2!=nullptr){
        if(knearest_2->Kmagnitude-knearest_1->Kmagnitude>KED_KMAG_LIMIT&&vector_angle(knearest_2->K, knearest_1->K)<=HALF_PI&&vector_angle(knearest_2->K, knearest_1->K)>1.0e-6) break;
        knearest_2=knearest_2->next;
    }
}

void KED::rotate_by_first_knearest(int zone[3])
{
    axes[2][0]=double(zone[0]); axes[2][1]=double(zone[1]); axes[2][2]=double(zone[2]);
    vector_normalize(axes[2], axes[2]);
    vector_zero(axes[2], axes[2]);
    vector_copy(axes[0], knearest_1->K);
    vector_normalize(axes[0], axes[0]);
    vector_zero(axes[0], axes[0]);
    vector_cross(axes[1], axes[0], axes[2]);
    vector_normalize(axes[1], axes[1]);
    vector_zero(axes[1], axes[1]);
    vector_cross(axes[0], axes[1], axes[2]);
    vector_normalize(axes[0], axes[0]);
    vector_zero(axes[0], axes[0]);

}

void KED::rotate(double x[3], double y[3])
{
    vector_copy(axes[0], x);
    vector_normalize(axes[0], axes[0]);
    vector_zero(axes[0], axes[0]);
    vector_copy(axes[1], y);
    vector_normalize(axes[1], axes[1]);
    vector_zero(axes[1], axes[1]);
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%d %d %d], y-[%d %d %d]", x[0], x[1], x[2], y[0], y[1], y[2]);
        exit(1);
    }
}

void KED::ked(char *ked_path)
{
    FILE *fp=nullptr;
    fp=fopen(ked_path,"w");
    fprintf(fp, "# N_1\tN_2\tN_3\tK_1\tK_2\tK_3\tx\ty\tz\tintensity\tintensity_norm (%d points, rotated by x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], and z-[%.8f %.8f %.8f])\n", numk-1,
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double *pos_x=nullptr, *pos_y=nullptr, *intensity=nullptr;
    callocate(&pos_x, numk-1, 0.0); 
    callocate(&pos_y, numk-1, 0.0); 
    callocate(&intensity, numk-1, 0.0);
    KED_KNODE *ktemp=khead->next;
    double constn=100.0/intensity_max;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        double intensity_norm=constn*ktemp->intensity;
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), ktemp->K[0], ktemp->K[1], ktemp->K[2], xyz[0], xyz[1], xyz[2], ktemp->intensity, intensity_norm);
        fflush(fp);
        pos_x[i-1]=xyz[0]; pos_y[i-1]=xyz[1]; intensity[i-1]=intensity_norm;
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", ked_path);

    char png_path[strlen(ked_path)+5];
    strcpy(png_path, ked_path); strcat(png_path, ".png");
    img(png_path, pos_x, pos_y, intensity, numk-1, Kmagnitude_max);
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void KED::ked(char *ked_path, double sigma, double dx)
{
    int    nbin=2*round(Kmagnitude_max/dx)+1;
    int    nbin_half=nbin/2;
    double *pos_x=nullptr, **intensity=nullptr;
    callocate(&pos_x, nbin, 0.0);
    callocate_2d(&intensity, nbin, nbin, 0.0);
    for(int i=0;i<=nbin_half;i++){
        pos_x[nbin_half+i]=double(i)*dx; 
        pos_x[nbin_half-i]=-double(i)*dx;
    }

    double **intensity_c=nullptr; 
    callocate_2d(&intensity_c, nbin, nbin, 0.0);
    KED_KNODE *ktemp=khead->next;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        int m=round(xyz[0]/dx)+nbin_half;
        int n=round(xyz[1]/dx)+nbin_half;
        gaussian(intensity_c, pos_x, pos_x, nbin, ktemp->intensity, pos_x[m], pos_x[n], sigma);
        for(int j=0;j<nbin;j++){
            for(int k=0;k<nbin;k++){
                intensity[j][k]+=intensity_c[j][k];
                intensity_c[j][k]=0.0;
            }
        }
        ktemp=ktemp->next;
    }
    deallocate_2d(intensity_c, nbin);

    double imax=0.0, imin=1.0e8;
    for(int i=0;i<nbin;i++){
        for(int j=0;j<nbin;j++){
            if(imax<intensity[i][j]) imax=intensity[i][j];
            if(imin>intensity[i][j]) imin=intensity[i][j];
        }
    }
    printf("[INFO] Number of profiled diffraction intensity: %d\n", nbin*nbin);
    printf("[INFO] Range of profiled diffraction intensity: %.8f %.8f\n", imin, imax);

    FILE *fp=nullptr;
    fp=fopen(ked_path,"w");
    fprintf(fp, "# x\ty\tintensity\tintensity_norm (%d points, rotated by x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], and z-[%.8f %.8f %.8f])\n", nbin*nbin,
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double constn=100.0/imax;
    for(int i=0;i<nbin;i++){
        for(int j=0;j<nbin;j++){
            fprintf(fp, "%.8f\t%.8f\t%.8f\t", pos_x[j], pos_x[i], intensity[i][j]);
            intensity[i][j]*=constn;
            fprintf(fp, "%.8f\n", intensity[i][j]);
            fflush(fp);
        }
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", ked_path);

    char png_path[strlen(ked_path)+5];
    strcpy(png_path, ked_path); strcat(png_path, ".png");
    img(png_path, intensity, nbin);
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void KED::img(char *png_path, double *x, double *y, double *value, int num, double limit)
{
    double height=6.0, width=6.0;
    int tick_max=int(limit);
    int n_major_tick=2*tick_max+1;
    double *major_ticks; mallocate(&major_ticks, n_major_tick);
    for(int i=0;i<n_major_tick;i++){
        major_ticks[i]=double(-tick_max+i);
    }
    GRAPH graph(width, height, 300);
    graph.set_xlim(-limit, limit);
    graph.set_ylim(-limit, limit);
    graph.set_xticks(major_ticks, n_major_tick);
    graph.set_yticks(major_ticks, n_major_tick);
    graph.set_tick_in(false);
    graph.scatter(x, y, value, num);
    graph.draw(png_path);
}

void KED::img(char *png_path, double **value, int nump)
{
    image_array(png_path, value, nump, nump, 6.0, 6.0, 300, false);
}

KED::KED(char *ked3_path)
{
    FILE *fp=fopen(ked3_path, "r");
    if(fp==NULL){
        printf("[ERROR] Unable to open file %s\n", ked3_path);
        exit(1);
    }
    fseek(fp, 0, SEEK_SET);

    char readbuff[MAX_LENTH_IN_NAME];
    int  keyi=0;
    for(int i=0;(keyi<2)&&i<MAX_WORD_NUMBER;i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, "ELECTRON_WAVELENTH")){
            fscanf(fp, "%lf", &lambda);
            keyi++;
        }
        if(0==strcmp(readbuff, "DIFFRACTION_INTENSITY")){
            int num;
            fscanf(fp, "%d", &num);
            double hkl[3], K[3];
            double Kmagnitude, intensity;
            for(int i=0;i<num;i++){
                fscanf(fp, "%lf", &hkl[0]);
                fscanf(fp, "%lf", &hkl[1]);
                fscanf(fp, "%lf", &hkl[2]);
                fscanf(fp, "%lf", &K[0]);
                fscanf(fp, "%lf", &K[1]);
                fscanf(fp, "%lf", &K[2]);
                fscanf(fp, "%lf", &Kmagnitude);
                fscanf(fp, "%lf", &intensity);
                add_k_node(hkl, K, Kmagnitude, intensity);
                if(Kmagnitude_max<Kmagnitude) Kmagnitude_max=Kmagnitude;
            }
            keyi++;
        }
    }
    if(keyi<2){
        printf("[ERROR] Unable to read file %s\n", ked3_path);
        exit(1);
    }
}

void KED::ked3(char *ked3_path)
{
    FILE *fp=fopen(ked3_path, "w");
    fprintf(fp, "ELECTRON_WAVELENTH\n");
    fprintf(fp, "%.8f\n", lambda);
    fprintf(fp, "DIFFRACTION_INTENSITY\n");
    fprintf(fp, "%d\n", numk-1);
    KED_KNODE *ktemp=khead->next;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->Kmagnitude, ktemp->intensity);
        fflush(fp);
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Three-dimensional %d diffraction intensity for calculating Kikuchi pattern stored in %s\n", numk-1, ked3_path);
}

// void KED::vtk(char *vtk_path)
// {
//     int nump, dimension[3];
//     for(int i=0;i<3;i++){
//         dimension[i]=kmax[i]-kmin[i]+1;
//     }
//     FILE *fp=nullptr;
//     fp=fopen(vtk_path,"w");
//     fprintf(fp, "# vtk DataFile Version 3.0\n");
//     fprintf(fp, "ELECTRON DIFFRACTION INTENSITY DISTRUBUTION\n");
//     fprintf(fp, "ASCII\n");
//     fprintf(fp, "DATASET STRUCTURED_POINTS\n");
//     fprintf(fp, "DIMENSIONS %d %d %d\n", dimension[0],  dimension[1], dimension[2]);
//     fprintf(fp, "ASPECT_RATIO %g %g %g\n", spacingK[0], spacingK[1], spacingK[2]);
//     fprintf(fp, "ORIGIN %g %g %g\n", kmin[0]*spacingK[0], kmin[1]*spacingK[1], kmin[2]*spacingK[2]);
//     fprintf(fp, "POINT_DATA %d\n", dimension[0]*dimension[1]*dimension[2]);
//     fprintf(fp, "SCALARS intensity float\n");
//     fprintf(fp, "LOOKUP_TABLE default\n");
//     double ***data; 
//     callocate_3d(&data, dimension[0], dimension[1], dimension[2], -1.0);
//     KED_KNODE *ktemp=khead;
//     for(int i=0;i<numk&&ktemp!=nullptr;i++){
//         int ix=ktemp->hkl[0]-kmin[0];
//         int iy=ktemp->hkl[1]-kmin[1];
//         int iz=ktemp->hkl[2]-kmin[2];
//         data[ix][iy][iz]=ktemp->intensity;
//         ktemp=ktemp->next;
//     }
//     for(int il=0;il<dimension[2];il++){
//         for(int ik=0;ik<dimension[1];ik++){
//             for(int ih=0;ih<dimension[0];ih++){
//                 fprintf(fp, "%g\n", data[ih][ik][il]);
//                 fflush(fp);
//             }
//         }
//     }
//     deallocate_3d(data, dimension[0], dimension[1]);
//     fclose(fp);
//     printf("[INFO] Visualized data for three-dimensional kinematical electron pattern stored in %s.\n", vtk_path);
// }

KKD::KKD(char *ked3_path, double xaxis[3], double yaxis[3], double zaxis[3], double thickness, double threshold, double ratiox, double ratioy, int npx, int npy, bool is_stereo_proj)
{
    printf("[INFO] Starting computation of kinematical Kikuchi diffraction...\n");
    KED ked(ked3_path);
    double kn=1.0/ked.lambda;
    printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", ked.lambda);
    printf("[INFO] Number of diffraction intensity: %d\n", ked.numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", ked.intensity_min, ked.intensity_max);
    ked.filter_diffraction_intensity(threshold);
    printf("[INFO] Number of filtered diffraction intensity: %d\n", ked.numk);
    printf("[INFO] Range of filtered diffraction intensity: %.8f %.8f\n", ked.intensity_min, ked.intensity_max);
    compute_Kikuchi_sphere_projection(xaxis, yaxis, zaxis, kn, ratiox, ratioy, npx, npy, is_stereo_proj);
    printf("[INFO] Kikuchi pattern has %d pixels along x-[%.8f %.8f %.8f] and %d pixels along y-[%.8f %.8f %.8f] under zone-[%.8f %.8f %.8f]\n", 
            numpx, axes[0][0], axes[0][1], axes[0][2], numpy, axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    printf("[INFO] Kikuchi pattern has %.8f distance [degree] along x axis and %.8f distance [degree] along y axis\n", thetax*RAD_TO_DEG, thetay*RAD_TO_DEG);
    printf("[INFO] Kikuchi pattern has %.8f distance [Angstrom-1] along x axis and %.8f distance [Angstrom-1] along y axis\n", thetax*kn, thetay*kn);
    compute_Kikuchi_intensity_projection(&ked, thickness, kn);
    printf("[INFO] Number of Kikuchi band on Kikuchi pattern: %d\n", numk);
    printf("[INFO] Range of Kikuchi intensity on Kikuchi pattern: %.8f %.8f\n", intensity_min, intensity_max);
    KKD_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        printf("[INFO] Kikuchi band %d: N-[%d %d %d] K-[%.8f %.8f %.8f] Kwidth-%.8f Kintensity1-%.8f Kintensity2-%.8f\n", i+1, 
        int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->Kwidth, ktemp->intensity1, ktemp->intensity2);
        ktemp=ktemp->next;
    }
    printf("[INFO] Ending computation of kinematical Kikuchi diffraction\n");
}

KKD::~KKD()
{
    if(0!=numpx&&0!=numpy){
        deallocate_3d(screenK0, numpy, numpx);
        deallocate_2d(screenI, numpy);
    }
    KKD_KNODE *cur=khead;
    while(cur!=nullptr){
        KKD_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void KKD::compute_Kikuchi_sphere_projection(double xaxis[3], double yaxis[3], double zaxis[3], double kn, double ratiox, double ratioy, int npx, int npy, bool is_stereo_proj)
{
    vector_copy(axes[0], xaxis);
    vector_normalize(axes[0], axes[0]);
    vector_copy(axes[1], yaxis);
    vector_normalize(axes[1], axes[1]);
    vector_copy(axes[2], zaxis);
    vector_normalize(axes[2], axes[2]);
    for(int i=0;i<3;i++){
        vector_zero(axes[i], axes[i]);
    }
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], z-[%.8f %.8f %.8f]", 
                xaxis[0], xaxis[1], xaxis[2], yaxis[0], yaxis[1], yaxis[2], zaxis[0], zaxis[1], zaxis[2]);
        exit(1);
    }

    numpx=npx; numpy=npy;
    if(0==numpx%2) numpx++;
    if(0==numpy%2) numpy++;
    callocate_3d(&screenK0, numpy, numpx, 3, 0.0);
    int impx=numpx/2, impy=numpy/2;
    int err;
    if(is_stereo_proj){
        for(int i=0;i<numpy;i++){
            for(int j=0;j<numpx;j++){
                double xy[2]={double(i-impy)/double(impy)*ratioy, double(j-impx)/double(impx)*ratiox};
                double xyz[3];
                compute_sphere_from_stereographic_projection(xyz, err, xy);
                if(0==err){
                    vector_transform(xyz, xyz, axes);
                    vector_copy(screenK0[i][j], xyz);
                }
            }
        }
    }else{
        for(int i=0;i<numpy;i++){
            for(int j=0;j<numpx;j++){
                double xy[2]={double(i-impy)/double(impy)*ratioy, double(j-impx)/double(impx)*ratiox};
                double xyz[3];
                compute_sphere_from_orthographic_projection(xyz, err, xy);
                if(0==err){
                    vector_transform(xyz, xyz, axes);
                    vector_copy(screenK0[i][j], xyz);
                }
            }
        }
    }
    thetax=acos(vector_dot(screenK0[impy][0], screenK0[impy][numpx-1]));
    thetay=acos(vector_dot(screenK0[0][impx], screenK0[numpy-1][impx]));
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            vector_constant(screenK0[i][j], kn, screenK0[i][j]);
        }
    }
}

void KKD::add_k_node(double hkl[3], double K[3], double Kwidth, double intensity1, double intensity2)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new KKD_KNODE;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kwidth=Kwidth;
        ktail->intensity1=intensity1;
        ktail->intensity2=intensity2;
        if(intensity_max<intensity1) intensity_max=intensity1;
        if(intensity_min>intensity1) intensity_min=intensity1;
        if(intensity_max<intensity2) intensity_max=intensity2;
        if(intensity_min>intensity2) intensity_min=intensity2;
        numk++;
    }else{
        ktail->next=new KKD_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kwidth=Kwidth;
        ktail->intensity1=intensity1;
        ktail->intensity2=intensity2;
        if(intensity_max<intensity1) intensity_max=intensity1;
        if(intensity_min>intensity1) intensity_min=intensity1;
        if(intensity_max<intensity2) intensity_max=intensity2;
        if(intensity_min>intensity2) intensity_min=intensity2;
        numk++;
    }
}

void KKD::compute_Kikuchi_intensity_projection(KED *ked, double thickness, double kn)
{
    callocate_2d(&screenI, numpy, numpx, 0.0);
    printf("[INFO] Starting projection of diffraction intensity on the Kikuchi pattern...\n");
    clock_t start, finish;
    start=clock();
    quick_sort(ked->khead);
    KED_KNODE *ktemp=ked->khead;
    int count=0;
    bool is_count=false;
    for(int i=0;i<ked->numk&&ktemp!=nullptr;i++){
        // double upper_bound=sqrt(kn*kn+ktemp->Kmagnitude*thickness/2.0);
        // double lower_bound=sqrt(kn*kn-ktemp->Kmagnitude*thickness/2.0);
        double upper_bound=ktemp->Kmagnitude+thickness/2.0;
        double lower_bound=ktemp->Kmagnitude-thickness/2.0;
        double hkl[3]; vector_copy(hkl, ktemp->hkl);
        double intensity=ktemp->intensity;
        for(int j=0;j<numpy;j++){
            for(int k=0;k<numpx;k++){
                // double d[3];
                // vector_difference(d, screenK0[j][k], ktemp->K);
                // double proj=vector_length(d);
                double proj=vector_dot(screenK0[j][k], ktemp->K)/ktemp->Kmagnitude;
                if(proj<=upper_bound&&proj>=lower_bound&&intensity>screenI[j][k]){
                    screenI[j][k]=intensity;
                    is_count=true;
                }
            }
        }
        ktemp=ktemp->next;
        if(is_count){
            if((0==int(hkl[0])+int(ktemp->hkl[0]))&&(0==int(hkl[1])+int(ktemp->hkl[1]))&&(0==int(hkl[2])+int(ktemp->hkl[2]))){
                for(int j=0;j<numpy;j++){
                    for(int k=0;k<numpx;k++){
                        // double d[3];
                        // vector_difference(d, screenK0[j][k], ktemp->K);
                        // double proj=vector_length(d);
                        double proj=vector_dot(screenK0[j][k], ktemp->K)/ktemp->Kmagnitude;
                        if(proj<=upper_bound&&proj>=lower_bound&&ktemp->intensity>screenI[j][k]){
                            screenI[j][k]=ktemp->intensity;
                        }
                    }
                }
                add_k_node(ktemp->hkl, ktemp->K, ktemp->Kmagnitude*2.0, ktemp->intensity, intensity);
            }
        }
        ktemp=ktemp->next;
        is_count=false;
        count++;
        if(0==count%1000){
            printf("[INFO] Completed diffraction intensity %d of %d\n", count, ked->numk);
        }
    }
    printf("[INFO] Ending projection of diffraction intensity on the Kikuchi pattern\n");
    finish=clock();
    printf("[INFO] Projection time [s]: %.8f.\n", double(finish-start)/CLOCKS_PER_SEC);
}

void KKD::kkd(char* kkd_path, char background)
{
    FILE *fp=nullptr;
    fp=fopen(kkd_path,"w");
    fprintf(fp, "KIKUCHI_IMAGE_SIZE\n");
    fprintf(fp, "%d %d\n", numpx, numpy);
    fprintf(fp, "KIKUCHI_IMAGE_VALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenI[i][j]);
            fflush(fp);
        }
    }
    printf("[INFO] Information for Kikuchi pattern stored in %s\n", kkd_path);
    char png_path[strlen(kkd_path)+5];
    strcpy(png_path, kkd_path); strcat(png_path, ".png");
    img(png_path, intensity_max, 0.0, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
}

void KKD::kkd(char* kkd_path, double vmax, double vmin, char background)
{
    FILE *fp=nullptr;
    fp=fopen(kkd_path,"w");
    fprintf(fp, "KIKUCHI_IMAGE_SIZE\n");
    fprintf(fp, "%d %d\n", numpx, numpy);
    fprintf(fp, "KIKUCHI_IMAGE_VALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenI[i][j]);
            fflush(fp);
        }
    }
    printf("[INFO] Information for Kikuchi pattern stored in %s\n", kkd_path);
    char png_path[strlen(kkd_path)+5];
    strcpy(png_path, kkd_path); strcat(png_path, ".png");
    img(png_path, vmax, vmin, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
}

void KKD::img(char* png_path, double vmax, double vmin, char background)
{
    double *wdata;
    unreshape_2d(&wdata, screenI, numpy, numpx);
    unsigned char *pixels;
    int num=numpx*numpy;
    mallocate(&pixels, 3*num);

    double vdiff=vmax-vmin;
    unsigned char rgb_min[3]={0}, rgb_max[3]={255};
    switch(background)
    {
    case 'b':
        rgb_min[0]=rgb_min[1]=rgb_min[2]=0;
        rgb_max[0]=rgb_max[1]=rgb_max[2]=255;
        break;
    case 'w':
        rgb_min[0]=rgb_min[1]=rgb_min[2]=255;
        rgb_max[0]=rgb_max[1]=rgb_max[2]=0;
        break;
    default:
        printf("[ERROR] Unrecognized background %s\n", background);
        exit(1);
    }
    unsigned char rgb_diff[3]; vector_difference(rgb_diff, rgb_max, rgb_min);
    unsigned char rgb[3];
    for(int i=0;i<num;i++){
        if(wdata[i]>=vmax){
            vector_copy(rgb, rgb_max);
        }else if(wdata[i]<=vmin){
            vector_copy(rgb, rgb_min);
        }else{
            rgb[0]=rgb_min[0]+int((wdata[i]-vmin)/vdiff*rgb_diff[0]);
            rgb[1]=rgb_min[1]+int((wdata[i]-vmin)/vdiff*rgb_diff[1]);
            rgb[2]=rgb_min[2]+int((wdata[i]-vmin)/vdiff*rgb_diff[2]);
        }
        pixels[i*3]=rgb[0]; pixels[i*3+1]=rgb[1]; pixels[i*3+2]=rgb[2];
    }
    image_pixels(png_path, pixels, numpx, numpy);
}