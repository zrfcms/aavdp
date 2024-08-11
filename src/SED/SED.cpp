#include "SED.h"

SED::SED(EMODEL *model, double Kmagnitude_max, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of selected-area electron diffraction...\n");
    radiusK=Kmagnitude_max;
    radiusE=1.0/model->lambda;
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(radiusK/spacingK[i]);
    }
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    compute_diffraction_intensity(model);
    printf("[INFO] Ending computation of selected-area electron diffraction\n");
}

SED::SED(EMODEL *model, int zone[3], double thickness, double Kmagnitude_max, double spacing[3], bool is_spacing_auto)
{
    printf("[INFO] Starting computation of selected-area electron diffraction...\n");
    radiusK=Kmagnitude_max;
    radiusE=1.0/model->lambda;
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(radiusK/spacingK[i]);
    }
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    compute_diffraction_intensity(model, zone, thickness);
    printf("[INFO] Ending computation of selected-area electron diffraction\n");
}

SED::~SED()
{
    SED_KNODE *cur=khead;
    while(cur!=nullptr){
        SED_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

// void SED::count_diffraction_vectors(EMODEL *model)
// {
//     printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
//     printf("[INFO] Radius of reciprocal sphere (Angstrom-1): %.8f\n", radiusK);
//     int NspacingK[3];
//     for(int i=0;i<3;i++){
//         NspacingK[i]=ceil(radiusK/spacingK[i]);
//     }
//     vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
//     for(int ih=kmin[0];ih<=kmax[0];ih++){
//         for(int ik=kmin[1];ik<=kmax[1];ik++){
//             for(int il=kmin[2];il<=kmax[2];il++){
//                 double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
//                 double Kmag=model->get_reciprocal_vector_length(K);
//                 if(Kmag<radiusK){
//                     numk++;
//                 }
//             }
//         }
//     }
//     printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
//     printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
//     printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
//     printf("[INFO] Number of diffraction vector in reciprocal space: %d\n", numk);
// }

// void SED::count_diffraction_vectors(EMODEL *model, int zone[3], double thickness)
// {
//     printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
//     printf("[INFO] Radius of reciprocal sphere (Angstrom-1): %.8f\n", radiusK);
//     int NspacingK[3];
//     for(int i=0;i<3;i++){
//         NspacingK[i]=ceil(radiusK/spacingK[i]);
//     }
//     double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
//     vector_normalize(n_zone, n_zone);
//     vector_constant(n_zone, radiusE, n_zone);
//     double upper_bound=radiusE+thickness/2.0;
//     double lower_bound=radiusE-thickness/2.0;
//     printf("[INFO] Range along zone-[%d %d %d] in reciprocal space (Angstrom-1): %.8f %.8f\n", zone[0], zone[1], zone[2], lower_bound, upper_bound);
//     for(int ih=-NspacingK[0];ih<=NspacingK[0];ih++){
//         for(int ik=-NspacingK[1];ik<=NspacingK[1];ik++){
//             for(int il=-NspacingK[2];il<=NspacingK[2];il++){
//                 double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
//                 double Kmag=model->get_reciprocal_vector_length(K);
//                 if(Kmag<radiusK){
//                     double d[3];
//                     model->reciprocal_to_cartesian(K, K);
//                     vector_difference(d, K, n_zone);
//                     double dmag=vector_length(d);
//                     if((dmag>lower_bound)&&(dmag<upper_bound)){
//                         if(ih<kmin[0]) kmin[0]=ih;
//                         if(ik<kmin[1]) kmin[1]=ik;
//                         if(il<kmin[2]) kmin[2]=il;
//                         if(ih>kmax[0]) kmax[0]=ih;
//                         if(ik>kmax[1]) kmax[1]=ik;
//                         if(il>kmax[2]) kmax[2]=il;
//                         numk++;
//                     }
//                 }
//             }
//         }
//     }
//     printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
//     printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
//     printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
//     printf("[INFO] Number of diffraction vector in reciprocal space: %d\n", numk);
// }

void SED::add_k_node(int h, int k, int l, double K[3], double Kmag, double intensity)
{
    if(ktail==nullptr){
        khead=ktail=new SED_KNODE;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        vector_copy(ktail->K, K);
        ktail->Kmag=Kmag;
        ktail->intensity=intensity;
        numk++;
    }else{
        ktail->next=new SED_KNODE;
        ktail=ktail->next;
        ktail->k[0]=h; ktail->k[1]=k; ktail->k[2]=l;
        vector_copy(ktail->K, K);
        ktail->Kmag=Kmag;
        ktail->intensity=intensity;
        numk++;
    }
}

void SED::compute_diffraction_intensity(EMODEL *model)
{
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int num=(kmax[0]-kmin[0]+1)*(kmax[1]-kmin[1]+1)*(kmax[2]-kmin[2]+1);
    int count=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(Kmag<radiusK){
                    model->reciprocal_to_cartesian(K, K);
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, K);
                    if(intensity>ZERO_LIMIT){
                        if(Kmag<ZERO_LIMIT){
                            intensity_0=intensity;
                        }else{
                            if(intensity>intensity_max) intensity_max=intensity;
                            if(intensity<intensity_min) intensity_min=intensity;
                            add_k_node(ih, ik, il, K, Kmag, intensity);
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
    printf("[INFO] Computation time [s]: %.2f\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Intensity at the transmission spot: %.8f\n", intensity_0);
    printf("[INFO] Number of diffraction intensity: %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void SED::compute_diffraction_intensity(EMODEL *model, int zone[3], double thickness)
{
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    vector_constant(n_zone, radiusE, n_zone);
    double upper_bound=radiusE+thickness/2.0;
    double lower_bound=radiusE-thickness/2.0;
    printf("[INFO] Range along zone-[%d %d %d] in reciprocal space (Angstrom-1): %.8f %.8f\n", zone[0], zone[1], zone[2], lower_bound, upper_bound);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int num=(kmax[0]-kmin[0]+1)*(kmax[1]-kmin[1]+1)*(kmax[2]-kmin[2]+1);
    int count=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(Kmag<radiusK){
                    double d[3];
                    model->reciprocal_to_cartesian(K, K);
                    vector_difference(d, K, n_zone);
                    double dmag=vector_length(d);
                    if((dmag>lower_bound)&&(dmag<upper_bound)){
                        double theta=asin(0.5*model->lambda*Kmag);
                        double intensity=model->get_diffraction_intensity(theta, K);
                        if(intensity>ZERO_LIMIT){
                            if(Kmag<ZERO_LIMIT){
                                intensity_0=intensity;
                            }else{
                                if(intensity>intensity_max) intensity_max=intensity;
                                if(intensity<intensity_min) intensity_min=intensity;
                                add_k_node(ih, ik, il, K, Kmag, intensity);
                            }
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
    printf("[INFO] Intensity at the transmission spot: %.8f\n", intensity_0);
    printf("[INFO] Number of diffraction intensity: %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void SED::vtk(const char *vtk_path, double threshold)
{
    int nump, dimension[3];
    for(int i=0;i<3;i++){
        dimension[i]=kmax[i]-kmin[i]+1;
    }
    FILE *fp=nullptr;
    fp=fopen(vtk_path,"w");
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "ELECTRON DIFFRACTION INTENSITY DISTRUBUTION\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", dimension[0],  dimension[1], dimension[2]);
    fprintf(fp, "ASPECT_RATIO %g %g %g\n", spacingK[0], spacingK[1], spacingK[2]);
    fprintf(fp, "ORIGIN %g %g %g\n", kmin[0]*spacingK[0], kmin[1]*spacingK[1], kmin[2]*spacingK[2]);
    fprintf(fp, "POINT_DATA %d\n", dimension[0]*dimension[1]*dimension[2]);
    fprintf(fp, "SCALARS intensity float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    double ***data; 
    callocate_3d(&data, dimension[0], dimension[1], dimension[2], -1.0);
    SED_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        int ix=ktemp->k[0]-kmin[0];
        int iy=ktemp->k[1]-kmin[1];
        int iz=ktemp->k[2]-kmin[2];
        data[ix][iy][iz]=ktemp->intensity;
        ktemp=ktemp->next;
    }
    data[-kmin[0]][-kmin[1]][-kmin[2]]=intensity_0;
    for(int il=0;il<dimension[2];il++){
        for(int ik=0;ik<dimension[1];ik++){
            for(int ih=0;ih<dimension[0];ih++){
                if(data[ih][ik][il]<=threshold) data[ih][ik][il]=0.0;
                fprintf(fp, "%g\n", data[ih][ik][il]);
                fflush(fp);
            }
        }
    }
    deallocate_3d(data, dimension[0], dimension[1]);
    fclose(fp);
    printf("[INFO] Visualized data for three-dimensional selected-area electron pattern stored in %s.\n", vtk_path);
}

void SED::sed(const char *sed_path, double threshold)
{
    FILE *fp=nullptr;
    fp=fopen(sed_path,"w");
    int numi=0;
    SED_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        if(ktemp->intensity>threshold) numi++;
        ktemp=ktemp->next;
    }

    fprintf(fp, "# x\ty\tz\tintensity\tintensity_norm (%d points)\n", numi+1);
    fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", 0.0, 0.0, 0.0, intensity_0, 100.0);
    double consti=100.0/intensity_max;
    ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        if(ktemp->intensity>threshold){
            fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->intensity, consti*ktemp->intensity);
            fflush(fp);
        }
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for three-dimensional selected-area electron pattern stored in %s\n", sed_path);
}

void SED::sed(const char *sed_path, char *png_path, int xaxis[3], int yaxis[3], int zaxis[3], double threshold)
{
    double rotation[3][3];
    double z[3]={double(zaxis[0]), double(zaxis[1]), double(zaxis[2])};
    vector_normalize(z, z); vector_copy(rotation[2], z);
    double x[3]={double(xaxis[0]), double(xaxis[1]), double(xaxis[2])};
    if(0==xaxis[0]&&0==xaxis[1]&&0==xaxis[2]){
        double Kmin=1.0e8;
        SED_KNODE *ktemp=khead;
        for(int i=0;i<numk&&ktemp!=nullptr;i++){
            if(ktemp->Kmag<Kmin&&ktemp->intensity>threshold){
                vector_copy(x, ktemp->K);
                Kmin=ktemp->Kmag;
            }
            ktemp=ktemp->next;
        }
    }
    vector_normalize(x, x); vector_copy(rotation[0], x);
    double y[3]={double(yaxis[0]), double(yaxis[1]), double(yaxis[2])};
    if(0==yaxis[0]&&0==yaxis[1]&&0==yaxis[2]){
        vector_cross(y, x, z);
    }
    vector_normalize(y, y); vector_copy(rotation[1], y);

    FILE *fp=nullptr;
    fp=fopen(sed_path,"w");
    int numi=0;
    SED_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        if(ktemp->intensity>threshold) numi++;
        ktemp=ktemp->next;
    }
    fprintf(fp, "# x\ty\tz\tintensity\tintensity_norm (%d points, x-[%.2f %.2f %.2f], y-[%.2f %.2f %.2f], z-[%.2f %.2f %.2f])\n", numi+1, x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]);
    double *kvector_x, *kvector_y, *kintensity;
    callocate(&kvector_x, numi+1, 0.0); 
    callocate(&kvector_y, numi+1, 0.0); 
    callocate(&kintensity, numi+1, 0.0);
    int counti=0;
    kvector_x[counti]=kvector_y[counti]=0.0; kintensity[counti]=100.0;
    fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", 0.0, 0.0, 0.0, ktemp->intensity, 100.0);
    counti++;
    double consti=100.0/intensity_max;
    ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        if(ktemp->intensity>threshold){
            double xyz[3], intensity;
            vector_rotate(xyz, rotation, ktemp->K);
            intensity=consti*ktemp->intensity;
            kvector_x[counti]=xyz[0]; kvector_y[counti]=xyz[1]; kintensity[counti]=intensity;
            fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", xyz[0], xyz[1], xyz[2], ktemp->intensity, intensity);
            counti++;
            fflush(fp);
        }
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for selected-area electron pattern stored in %s.\n", sed_path);

    double height=6.0, width=6.0;
    int tick_max=int(floor(radiusK));
    int n_major_tick=2*tick_max+1;
    double *major_ticks; mallocate(&major_ticks, n_major_tick);
    for(int i=0;i<n_major_tick;i++){
        major_ticks[i]=double(-tick_max+i);
    }
    GRAPH graph(width, height, 300);
    graph.set_xlim(-radiusK, radiusK);
    graph.set_ylim(-radiusK, radiusK);
    graph.set_xticks(major_ticks, n_major_tick);
    graph.set_yticks(major_ticks, n_major_tick);
    graph.set_tick_in(false);
    graph.scatter(kvector_x, kvector_y, kintensity, numi);
    graph.draw(png_path);
    printf("[INFO] Image for selected-area electron pattern stored in %s\n", png_path);
}

void SED::restart(const char *restart_path)
{
    FILE *fp=fopen(restart_path, "w");
    fprintf(fp, "Kikuchi_sphere_radius\n");
    fprintf(fp, "%.8f\n", radiusE);
    fprintf(fp, "Kikuchi_line_intensity\n");
    fprintf(fp, "%d\n", numk);
    SED_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\n", ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->intensity);
        fflush(fp);
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Number of Kikuchi lines stored: %d\n", numk);
    printf("[INFO] Range of Kikuchi line intensity stored: %.8f %.8f\n", intensity_min, intensity_max);
    printf("[INFO] Restart data for Kikuchi pattern stored in %s\n", restart_path);
}