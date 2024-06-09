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
    count_diffraction_vectors(model);
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
    count_diffraction_vectors(model, zone, thickness);
    compute_diffraction_intensity(model, zone, thickness);
    printf("[INFO] Ending computation of selected-area electron diffraction\n");
}

SED::SED(const char *hdf5_path)
{
    size_t size1, size2;
    double *spacing, *Irange;
    int    **krange;
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.read("/ElectronDiffraction/ReciprocalSphereRadius", radiusK);
    hdf.read("/ElectronDiffraction/EwaldSphereRadius", radiusE);
    hdf.read_array("/ElectronDiffraction/ReciprocalSpacing", &spacing, size1);
    hdf.read_array("/ElectronDiffraction/DiffractionIntensityRange", &Irange, size1);
    hdf.read_array_2d("/ElectronDiffraction/DiffractionIndiceRange", &krange, size1, size2);
    hdf.read("/ElectronDiffraction/DiffractionVectorNumber", numk);
    hdf.read_array_2d("/ElectronDiffraction/DiffractionIndiceList", &kvectors, size1, size2);
    hdf.read_array_2d("/ElectronDiffraction/DiffractionVectorList", &Kvectors, size1, size2);
    hdf.read_array("/ElectronDiffraction/DiffractionIntensityList", &Kintensity, size1);
    hdf.close();
    for(int i=0;i<3;i++){
        spacingK[i]=spacing[i];
        kmin[i]=krange[i][0]; kmax[i]=krange[i][1];
    }
    intensity_min=Irange[0]; intensity_max=Irange[1];
    deallocate(spacing); deallocate(Irange); deallocate_2d(krange, 3);
}

SED::~SED()
{
    if(0!=numk){
        deallocate_2d(kvectors, numk);
        deallocate_2d(Kvectors, numk);
        deallocate(Kintensity);
    }
}

void SED::count_diffraction_vectors(EMODEL *model)
{
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Radius of reciprocal sphere (Angstrom-1): %.8f\n", radiusK);
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(radiusK/spacingK[i]);
    }
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(Kmag<radiusK){
                    numk++;
                }
            }
        }
    }
    printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("[INFO] Number of diffraction vector in reciprocal space: %d\n", numk);
}

void SED::count_diffraction_vectors(EMODEL *model, int zone[3], double thickness)
{
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
    printf("[INFO] Radius of reciprocal sphere (Angstrom-1): %.8f\n", radiusK);
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(radiusK/spacingK[i]);
    }
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    vector_constant(n_zone, radiusE, n_zone);
    double upper_bound=radiusE+thickness/2.0;
    double lower_bound=radiusE-thickness/2.0;
    printf("[INFO] Range along zone-[%d %d %d] in reciprocal space (Angstrom-1): %.8f %.8f\n", zone[0], zone[1], zone[2], lower_bound, upper_bound);
    for(int ih=-NspacingK[0];ih<=NspacingK[0];ih++){
        for(int ik=-NspacingK[1];ik<=NspacingK[1];ik++){
            for(int il=-NspacingK[2];il<=NspacingK[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(Kmag<radiusK){
                    double d[3];
                    model->reciprocal_to_cartesian(K, K);
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
    printf("[INFO] Range of Miller index h in reciprocal space: %d %d\n", kmin[0], kmax[0]);
    printf("[INFO] Range of Miller index k in reciprocal space: %d %d\n", kmin[1], kmax[1]);
    printf("[INFO] Range of Miller index l in reciprocal space: %d %d\n", kmin[2], kmax[2]);
    printf("[INFO] Number of diffraction vector in reciprocal space: %d\n", numk);
}

void SED::compute_diffraction_intensity(EMODEL *model)
{
    callocate_2d(&kvectors, numk, 3, 0);
    callocate_2d(&Kvectors, numk, 3, 0.0);
    callocate(&Kintensity, numk, 0.0);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmag=model->get_reciprocal_vector_length(K);
                if(Kmag<radiusK){
                    model->reciprocal_to_cartesian(K, K);
                    double theta=asin(0.5*model->lambda*Kmag);
                    double intensity=model->get_diffraction_intensity(theta, K);
                    if(intensity>intensity_max) intensity_max=intensity;
                    if(intensity<intensity_min) intensity_min=intensity;
                    kvectors[countk][0]=ih; kvectors[countk][1]=ik; kvectors[countk][2]=il;
                    Kvectors[countk][0]=K[0]; Kvectors[countk][1]=K[1]; Kvectors[countk][2]=K[2];
                    Kintensity[countk]=intensity;
                    countk++;
                    if(0==countk%1000){
                        printf("[INFO] Completed diffraction intensity %d of %d\n", countk, numk);
                    }
                }
            }
        }
    }
    printf("[INFO] Ending computation of diffraction intensity\n");
    finish=clock();
    printf("[INFO] Computation time [s]: %.2f\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Number of diffraction intensity: %d\n", countk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void SED::compute_diffraction_intensity(EMODEL *model, int zone[3], double thickness)
{
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    vector_constant(n_zone, radiusE, n_zone);
    double upper_bound=radiusE+thickness/2.0;
    double lower_bound=radiusE-thickness/2.0;
    callocate_2d(&kvectors, numk, 3, 0);
    callocate_2d(&Kvectors, numk, 3, 0.0);
    callocate(&Kintensity, numk, 0.0);
    printf("[INFO] Starting computation of diffraction intensity...\n");
    clock_t start, finish;
    start=clock();
    int countk=0;
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
                        if(intensity>intensity_max) intensity_max=intensity;
                        if(intensity<intensity_min) intensity_min=intensity;
                        kvectors[countk][0]=ih; kvectors[countk][1]=ik; kvectors[countk][2]=il;
                        Kvectors[countk][0]=K[0]; Kvectors[countk][1]=K[1]; Kvectors[countk][2]=K[2];
                        Kintensity[countk]=intensity;
                        countk++;
                        if(0==countk%1000){
                            printf("[INFO] Completed diffraction intensity %d of %d\n", countk, numk);
                        }
                    }
                }
            }
        }
    }
    printf("[INFO] Ending computation of diffraction intensity\n");
    finish=clock();
    printf("[INFO] Computation time [s]: %.2f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Number of diffraction intensity: %d\n", countk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
}

void SED::hdf5(const char *hdf5_path)
{
    double *spacing; mallocate(&spacing, 3);
    int    **krange; mallocate_2d(&krange, 3, 2);
    double *Irange; mallocate(&Irange, 2);
    for(int i=0;i<3;i++){
        spacing[i]=spacingK[i];
        krange[i][0]=kmin[i]; krange[i][1]=kmax[i];
        Irange[0]=intensity_min; Irange[1]=intensity_max;
    }
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.write_group("/ElectronDiffraction");
    hdf.write("/ElectronDiffraction/ReciprocalSphereRadius", radiusK);
    hdf.write("/ElectronDiffraction/EwaldSphereRadius", radiusE);
    hdf.write_array("/ElectronDiffraction/ReciprocalSpacing", spacing, 3);
    hdf.write_array_2d("/ElectronDiffraction/DiffractionIndiceRange", krange, 3, 2);
    hdf.write("/ElectronDiffraction/DiffractionVectorNumber", numk);
    hdf.write_array_2d("/ElectronDiffraction/DiffractionIndiceList", kvectors, numk, 3);
    hdf.write_array_2d("/ElectronDiffraction/DiffractionVectorList", Kvectors, numk, 3);
    hdf.write_array("/ElectronDiffraction/DiffractionIntensityList", Kintensity, numk);
    hdf.write_array("/ElectronDiffraction/DiffractionIntensityRange", Irange, 2);
    hdf.close();
    deallocate(spacing); deallocate(Irange); deallocate_2d(krange, 3);
    printf("[INFO] Information for selected-area electron pattern stored in %s.\n", hdf5_path);
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
                if(fabs(data[ih][ik][il])<1.0e-15) data[ih][ik][il]=0.0;
                fprintf(fp, "%g\n", data[ih][ik][il]);
                fflush(fp);
            }
        }
    }
    deallocate_3d(data, dimension[0], dimension[1]);
    fclose(fp);
    printf("[INFO] Visualized data for selected-area electron pattern stored in %s.\n", vtk_path);
}

void SED::sed(const char *sed_path, int xaxis[3], int zaxis[3], double thickness, double threshold_ratio)
{
    double rotation[3][3];
    double z[3]={double(zaxis[0]), double(zaxis[1]), double(zaxis[2])};
    vector_normalize(z, z); vector_copy(rotation[2], z);

    int numi=0; 
    int *Kindice; callocate(&Kindice, numk, 0);
    double ithreshold=threshold_ratio, imax=0.0;
    double imt=thickness/2.0;
    for(int i=0;i<numk;i++){
        if(Kintensity[i]>ithreshold&&fabs(vector_dot(Kvectors[i], z))<imt){
            Kindice[numi]=i;
            numi++;
            if(0==kvectors[i][0]&&0==kvectors[i][1]&&0==kvectors[i][2]) continue;
            if(imax<Kintensity[i]) imax=Kintensity[i];
        }
    }

    double x[3]={double(xaxis[0]), double(xaxis[1]), double(xaxis[2])};
    if(0==xaxis[0]&&0==xaxis[1]&&0==xaxis[2]){
        double Kmag_min=1.0e8;
        for(int i=0;i<numi;i++){
            if(0==kvectors[Kindice[i]][0]&&0==kvectors[Kindice[i]][1]&&0==kvectors[Kindice[i]][2]) continue;
            double Kmag=vector_length(Kvectors[Kindice[i]]);
            if(Kmag<Kmag_min){
                vector_copy(x, Kvectors[Kindice[i]]);
                Kmag_min=Kmag;
            }
        }
    }
    vector_normalize(x, x); vector_copy(rotation[0], x);

    double y[3];
    vector_cross(y, x, z);
    vector_normalize(y, y); vector_copy(rotation[1], y);

    FILE *fp=nullptr;
    fp=fopen(sed_path,"w");
    fprintf(fp, "# x\ty\tz\tintensity\tintensity_norm (%d points)\n", numi);
    double consti=100.0/imax;
    for(int i=0;i<numi;i++){
        if(0==kvectors[Kindice[i]][0]&&0==kvectors[Kindice[i]][1]&&0==kvectors[Kindice[i]][2]){
            fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", 0.0, 0.0, 0.0, Kintensity[Kindice[i]], 100.0);
            continue;
        }
        double xyz[3];
        vector_rotate(xyz, rotation, Kvectors[Kindice[i]]);
        fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", xyz[0], xyz[1], xyz[2], Kintensity[Kindice[i]], consti*Kintensity[Kindice[i]]);
        fflush(fp);
    }
    fclose(fp);
    printf("[INFO] Visualized data for selected-area electron pattern stored in %s.\n", sed_path);
}