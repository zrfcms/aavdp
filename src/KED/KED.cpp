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

KED::KED(MODEL *model, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of kinematical electron diffraction...\n");
    printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", model->lambda);
    Kmagnitude_max=Kmag_max;
    compute_diffraction_intensity(model, threshold, spacing, is_spacing_auto);
    filter_diffraction_intensity(threshold);
    printf("[INFO] Number of diffraction intensity: %d\n", numk);
    printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    printf("[INFO] Ending computation of kinematical electron diffraction\n");
}

KED::KED(MODEL *model, int zone[3], double thickness, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of kinematical electron diffraction...\n");
    printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", model->lambda);
    Kmagnitude_max=Kmag_max;
    compute_diffraction_intensity(model, zone, thickness, spacing, is_spacing_auto);
    filter_diffraction_intensity(threshold);
    printf("[INFO] Number of diffraction intensity: %d\n", numk);
    printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    find_first_and_second_knearests();
    if(knearest_2==nullptr){
        printf("[INFO] The first nearest diffraction vectors R1: [%.5f %.5f %.5f]\n", knearest_1->K[0], knearest_1->K[1], knearest_1->K[2]);
        printf("[WARN] Unable to find the second nearest diffraction vector\n");
    }else{
        printf("[INFO] The first and second nearest diffraction vectors R1, R2: [%.5f %.5f %.5f], [%.5f %.5f %.5f] (R2/R1 %.8f and angle %.8f)\n", 
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
        ktail->Kmagnitude=Kmagnitude;
        ktail->intensity=intensity;
        numk++;
    }else{
        ktail->next=new KED_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        ktail->Kmagnitude=Kmagnitude;
        ktail->intensity=intensity;
        if(intensity_max<intensity) intensity_max=intensity;
        if(intensity_min>intensity) intensity_min=intensity;
        numk++;
    }
}

void KED::compute_diffraction_intensity(MODEL *model, double threshold, double spacing[3], bool is_spacing_auto)
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
    printf("[INFO] Spacings along a*, b*, and c* in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Range along a*, b*, or c* in reciprocal space (Angstrom-1): %.8f %.8f\n", 0.0, Kmagnitude_max);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    double hkl[3];
    double K[3]={0.0};
    double intensity=model->get_diffraction_intensity(0.0, K, true);
    add_k_node(hkl, K, 0.0, intensity);
    clock_t start, finish;
    start=clock();
    int num=(2*kmax[0]+1)*(kmax[1]+1)*(2*kmax[2]+1);
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
                    double intensity=model->get_diffraction_intensity(Kmagnitude, K, false);
                    if(intensity>KED_INTENSITY_LIMIT){
                        add_k_node(hkl, K, Kmagnitude, intensity);
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
    printf("[INFO] Spacings along a*, b*, and c* in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Range along a*, b*, or c* in reciprocal space (Angstrom-1): %.8f %.8f\n", 0.0, Kmagnitude_max);
    printf("[INFO] Projection range along zone-[%.5f %.5f %.5f] in reciprocal space (Angstrom-1): %.8f %.8f\n", n_zone[0], n_zone[1], n_zone[2], lower_bound, upper_bound);
    double hkl[3]={0.0};
    double K[3]={0.0};
    double intensity=model->get_diffraction_intensity(0.0, K, true);
    add_k_node(hkl, K, 0.0, intensity);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int num=(kmax[0]-kmin[0]+1)*(kmax[1]-kmin[1]+1)*(kmax[2]-kmin[2]+1);
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
    vector_copy(axes[0], knearest_1->K);
    vector_normalize(axes[0], axes[0]);
    vector_cross(axes[1], axes[0], axes[2]);
    vector_normalize(axes[1], axes[1]);
    vector_cross(axes[0], axes[1], axes[2]);
    vector_normalize(axes[0], axes[0]);
}

void KED::rotate(int x[3], int y[3])
{
    axes[0][0]=double(x[0]); axes[0][1]=double(x[1]); axes[0][2]=double(x[2]);
    axes[1][0]=double(y[0]); axes[1][1]=double(y[1]); axes[1][2]=double(y[2]);
    vector_normalize(axes[0], axes[0]);
    vector_normalize(axes[1], axes[1]);
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%d %d %d], y-[%d %d %d]", x[0], x[1], x[2], y[0], y[1], y[2]);
        exit(1);
    }
}

void KED::ked(char *ked_path)
{
    FILE *fp=nullptr;
    fp=fopen(ked_path,"w");
    fprintf(fp, "# K_1\tK_2\tK_3\tx\ty\tz\tintensity\tintensity_norm (%d points, rotated by x-[%.5f %.5f %.5f], y-[%.5f %.5f %.5f], and z-[%.5f %.5f %.5f])\n", numk,
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double *pos_x=nullptr, *pos_y=nullptr, *intensity=nullptr;
    callocate(&pos_x, numk, 0.0); 
    callocate(&pos_y, numk, 0.0); 
    callocate(&intensity, numk, 0.0);
    KED_KNODE *ktemp=khead;
    fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", ktemp->K[0], ktemp->K[1], ktemp->K[2], 0.0, 0.0, 0.0, ktemp->intensity, 100.0);
    fflush(fp);
    pos_x[0]=pos_y[0]=0.0; intensity[0]=100.0;
    ktemp=ktemp->next;
    double constn=100.0/intensity_max;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        double intensity_norm=constn*ktemp->intensity;
        fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", ktemp->K[0], ktemp->K[1], ktemp->K[2], xyz[0], xyz[1], xyz[2], ktemp->intensity, intensity_norm);
        fflush(fp);
        pos_x[i]=xyz[0]; pos_y[i]=xyz[1]; intensity[i]=intensity_norm;
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", ked_path);

    char png_path[PATH_CHAR_NUMBER];
    strcpy(png_path, ked_path); strcat(png_path, ".png");
    img(png_path, pos_x, pos_y, intensity, numk, Kmagnitude_max);
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

// void KED::restart(const char *restart_path)
// {
//     FILE *fp=fopen(restart_path, "w");
//     fprintf(fp, "Kikuchi_sphere_radius\n");
//     fprintf(fp, "%.8f\n", radiusE);
//     fprintf(fp, "Kikuchi_line_intensity\n");
//     fprintf(fp, "%d\n", numk-1);
//     KED_KNODE *ktemp=khead->next;
//     for(int i=1;i<numk&&ktemp!=nullptr;i++){
//         fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\n", ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->intensity);
//         fflush(fp);
//         ktemp=ktemp->next;
//     }
//     fclose(fp);
//     printf("[INFO] Number of Kikuchi lines stored: %d\n", numk-1);
//     printf("[INFO] Range of Kikuchi line intensity stored: %.8f %.8f\n", intensity_min, intensity_max);
//     printf("[INFO] Restart data for Kikuchi pattern stored in %s\n", restart_path);
// }