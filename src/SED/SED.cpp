#include "SED.h"

SED::SED(EMODEL *model, double spacing[3], double Kmagnitude_max)
{
    printf("Starting...\n");
    spacingK[0]=spacing[0]; spacingK[1]=spacing[1]; spacingK[2]=spacing[2];
    count_diffraction_vectors(Kmagnitude_max);
    compute_diffraction_intensity(model, Kmagnitude_max);
    printf("Ending...\n");
}

SED::SED(EMODEL *model, double spacing[3], int zone[3], double thickness, double Kmagnitude_max)
{
    printf("Starting...\n");
    spacingK[0]=spacing[0]; spacingK[1]=spacing[1]; spacingK[2]=spacing[2];
    radiusE=1.0/model->lambda;
    count_diffraction_vectors(zone, thickness, Kmagnitude_max);
    compute_diffraction_intensity(model, zone, thickness, Kmagnitude_max);
    printf("Ending...\n");
}

SED::~SED()
{
    if(0!=numk){
        deallocate_2d(kvectors, numk);
        deallocate_2d(Kvectors, numk);
        deallocate(Kintensity);
    }
}

void SED::count_diffraction_vectors(double Kmagnitude_max)
{
    printf("Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    kmax[0]=NspacingK[0]; kmax[1]=NspacingK[1]; kmax[2]=NspacingK[2];
    kmin[0]=-NspacingK[0]; kmin[1]=-NspacingK[1]; kmin[2]=-NspacingK[2];
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3];
                K[0]=ih*spacingK[0];
                K[1]=ik*spacingK[1];
                K[2]=il*spacingK[2];
                double Kmag=vector_length(K);
                if(Kmag<Kmagnitude_max){
                    numk++;
                }
            }
        }
    }
    printf("Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("Number of diffraction vector in reciprocal space: %d\n", numk);
}

void SED::count_diffraction_vectors(int zone[3], double thickness, double Kmagnitude_max)
{
    printf("Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    n_zone[0]*=radiusE; n_zone[1]*=radiusE; n_zone[2]*=radiusE;
    double upper_bound=radiusE+thickness/2.0;
    double lower_bound=radiusE-thickness/2.0;
    printf("Range along zone-[%d %d %d] in reciprocal space (Angstrom-1): %.8f %.8f\n", lower_bound, upper_bound);
    for(int ih=-NspacingK[0];ih<=NspacingK[0];ih++){
        for(int ik=-NspacingK[1];ik<=NspacingK[1];ik++){
            for(int il=-NspacingK[2];il<=NspacingK[2];il++){
                double K[3];
                K[0]=ih*spacingK[0];
                K[1]=ik*spacingK[1];
                K[2]=il*spacingK[2];
                double Kmag=vector_length(K);
                if(Kmag<Kmagnitude_max){
                    double d[3];
                    vector_difference(d, K, n_zone);
                    double dmag=vector_length(d);
                    if((dmag>lower_bound)&&(dmag<upper_bound)){
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

void SED::compute_diffraction_intensity(EMODEL *model, double Kmagnitude_max)
{
    callocate_2d(&kvectors, numk, 3, 0);
    callocate_2d(&Kvectors, numk, 3, 0.0);
    callocate(&Kintensity, numk, 0.0);
    printf("Starting computation of diffraction intensity...\n");
    clock_t start, finish, time;
    start=clock();
    int countk=0;
    double intensity_min=1.0e8, intensity_max=0.0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3];
                K[0]=ih*spacingK[0];
                K[1]=ik*spacingK[1];
                K[2]=il*spacingK[2];
                double Kmag=vector_length(K);
                if(Kmag<Kmagnitude_max){
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, K);
                    if(intensity>intensity_max) intensity_max=intensity;
                    if(intensity<intensity_min) intensity_min=intensity;
                    kvectors[countk][0]=ih; kvectors[countk][1]=ik; kvectors[countk][2]=il;
                    Kvectors[countk][0]=K[0]; Kvectors[countk][1]=K[1]; Kvectors[countk][2]=K[2];
                    Kintensity[countk]=intensity;
                    countk++;
                    if(0==countk%1000){
                        printf("Completed diffraction intensity %d of %d\n", countk, numk);
                    }
                }
            }
        }
    }
    printf("Ending computation of diffraction intensity\n");
    finish=clock(); time=finish-start;
    printf("Computation time [s]: %.2f\n", time);
    printf("Number of diffraction intensity: %d\n", countk);
    printf("Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void SED::compute_diffraction_intensity(EMODEL *model, int zone[3], double thickness, double Kmagnitude_max)
{
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    n_zone[0]*=radiusE; n_zone[1]*=radiusE; n_zone[2]*=radiusE;
    double upper_bound=radiusE+thickness;
    double lower_bound=radiusE-thickness;
    callocate_2d(&kvectors, numk, 3, 0);
    callocate_2d(&Kvectors, numk, 3, 0.0);
    callocate(&Kintensity, numk, 0.0);
    printf("Starting computation of diffraction intensity...\n");
    clock_t start, finish, time;
    start=clock();
    int countk=0;
    double intensity_min=1.0e8, intensity_max=0.0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3];
                K[0]=ih*spacingK[0];
                K[1]=ik*spacingK[1];
                K[2]=il*spacingK[2];
                double Kmag=vector_length(K);
                if(Kmag<Kmagnitude_max){
                    double d[3];
                    vector_difference(d, K, n_zone);
                    double dmag=vector_length(d);
                    if((dmag>lower_bound)&&(dmag<upper_bound)){
                        double theta=asin(0.5*model->lambda*Kmag);
                        double intensity=model->get_diffraction_intensity(theta, K);
                        if(intensity>intensity_max) intensity_max=intensity;
                        if(intensity<intensity_min) intensity_min=intensity;
                        kvectors[countk][0]=ih; kvectors[countk][1]=ik; kvectors[countk][2]=il;
                        Kvectors[countk][0]=K[0]; Kvectors[countk][1]=K[1]; Kvectors[countk][2]=K[2];
                        Kintensity[countk]=intensity;
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
    finish=clock(); time=finish-start;
    printf("Computation time [s]: %.2f.\n", time);
    printf("Number of diffraction intensity: %d\n", countk);
    printf("Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void SED::vtk(const char *vtk_path)
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
    for(int i=0;i<numk;i++){
        int ix=kvectors[i][0]-kmin[0];
        int iy=kvectors[i][1]-kmin[1];
        int iz=kvectors[i][2]-kmin[2];
        data[ix][iy][iz]=Kintensity[i];
    }
    for(int il=0;il<dimension[2];il++){
        for(int ik=0;ik<dimension[1];ik++){
            for(int ih=0;ih<dimension[0];ih++){
                fprintf(fp, "%g\n", data[ih][ik][il]);
                fflush(fp);
            }
        }
    }
    fclose(fp);
    printf("Visualized data stored in %s.\n", vtk_path);
    deallocate_3d(data, dimension[0], dimension[1]);
}

void SED::kkd(const char *kkd_path)
{
    FILE *fp=nullptr;
    fp=fopen(kkd_path,"w");
    fprintf(fp,"# K_1\tK_2\tK_3\tIntensity\n");
    fprintf(fp,"Kikuchi_radius %.8f\n", radiusE);
    fprintf(fp,"Kvectors %d\n", numk);
    for(int i=0;i<numk;i++){
        fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\n", 
                    Kvectors[i][0], Kvectors[i][1], Kvectors[i][2], Kintensity[i]);
        fflush(fp);
    }
    fclose(fp);
    printf("Information for kinematic Kikuchi diffraction stored in %s.\n", kkd_path);
}