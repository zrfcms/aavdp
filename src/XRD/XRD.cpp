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

XRD::XRD(const char *xrd_path)
{
    FILE *fp=fopen(xrd_path, "r");
    if(fp==nullptr){
        printf("[ERROR] Unable to open file %s\n", xrd_path);
    }
    fseek(fp, 0, SEEK_SET);

    char readbuff[KEY_CHAR_NUMBER]; 
    int nkey=1, keyi=0;
    for(int i=0;(keyi<nkey)&&(i<KEY_LINE_NUMBER);i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, "xrd")){
            fscanf(fp, "%d", &numk);
            callocate_2d(&kvector, numk, 3, 0);
            callocate(&Ktheta, numk, 0.0);
            callocate(&Kfactor, numk, complex<double>(0.0, 0.0));
            callocate(&Kintensity, numk, 0.0);
            for(int j=0;j<numk;j++){
                double real=0.0, imag=0.0;
                fscanf(fp, "%d", &kvector[j][0]);
                fscanf(fp, "%d", &kvector[j][1]);
                fscanf(fp, "%d", &kvector[j][2]);
                fscanf(fp, "%lf", &Ktheta[j]);
                fscanf(fp, "%lf", &real); 
                fscanf(fp, "%lf", &imag);
                Kfactor[j].real(real); Kfactor[j].imag(imag);
                fscanf(fp, "%lf", &Kintensity[j]);
            }
            keyi++;
            continue;
        }
    }
    if(keyi!=nkey){
        printf("[ERROR] Unable to recognize information from %s\n", xrd_path);
    }
}

XRD::~XRD()
{
    if(0!=numk){
        deallocate_2d(kvector, numk);
        deallocate(Ktheta);
        deallocate(Kfactor);
        deallocate(Kintensity);
    }
}

void XRD::count_diffraction_vector(MODEL *model)
{
    int c_kmin[3]={kmin[0], kmin[1], kmin[2]}, c_kmax[3]={kmax[0], kmax[1], kmax[2]};
    for(int ih=c_kmin[0];ih<=c_kmax[0];ih++){
        for(int ik=c_kmin[1];ik<=c_kmax[1];ik++){
            for(int il=c_kmin[2];il<=c_kmax[2];il++){
                if(0==ih&&0==ik&&0==il) continue;
                double k[3]={double(ih), double(ik), double(il)};
                double Kmag=model->get_reciprocal_vector_length(k);
                if(2>=Kmag*model->lambda){
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
    printf("Number of diffraction vector in reciprocal space: %d\n", numk);
}

void XRD::compute_diffraction_intensity(MODEL *model)
{
    callocate_2d(&kvector, numk, 3, 0);
    callocate(&Ktheta, numk, 0.0);
    callocate(&Kfactor, numk, complex<double>(0.0, 0.0));
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
                    complex<double> factor=model->get_atomic_structure_factor(theta, k, is_lorentz_flag);
                    double intensity=model->get_diffraction_intensity(factor);
                    if(intensity<intensity_min) intensity_min=intensity;
                    if(intensity>intensity_max) intensity_max=intensity;
                    kvector[countk][0]=ih; kvector[countk][1]=ik; kvector[countk][2]=il;
                    Ktheta[countk]=theta; Kfactor[countk]=factor; Kintensity[countk]=intensity;
                    countk++;
                    if(0==countk%1000){
                        printf("Completed diffraction intensity %d of %d\n", countk, numk);
                    }
                }
            }
        }
    }
    printf("Ending computation of diffraction intensity\n");
    finish=clock();
    printf("Computation time [s]: %.8f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("Number of diffraction intensity: %d\n", countk);
    printf("Range of diffraction intensity: %.8f %.8f", intensity_min, intensity_max);
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
    callocate(&Kfactor, numk, complex<double>(0.0, 0.0));
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
                        Ktheta[countk]=theta; Kfactor[countk]=factor; Kintensity[countk]=intensity;
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
    printf("Range of diffraction intensity: %.8f %.8f", intensity_min, intensity_max);
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
    fprintf(fp,"# h\tk\tl\t2theta\tFreal\tFimag\tintensity\n");
    for(int i=0;i<numk;i++){
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\n", 
                    kvector[i][0], kvector[i][1], kvector[i][2], 2*Ktheta[i]/DEG_TO_RAD, 
                    Kfactor[i].real(), Kfactor[i].imag(), Kintensity[i]);
        fflush(fp);
    }
    fclose(fp);
    printf("Information for x-ray diffraction stored in %s.\n", xrd_path);
}

void XRD::unique()
{
    int    countk=0;
    int    **kuvector=nullptr;
    double *Kutheta=nullptr;
    complex<double> *Kufactor=nullptr;
    double *Kuintensity=nullptr;
    callocate_2d(&kuvector, numk, 3, 0);
    callocate(&Kutheta, numk, 0.0);
    callocate(&Kufactor, numk, complex<double>(0.0, 0.0));
    callocate(&Kuintensity, numk, 0.0);
    for(int i=0;i<numk;i++){
        if(i==0){
            vector_copy(kuvector[countk], kvector[i]);
            Kutheta[countk]=Ktheta[i]; Kufactor[countk]=Kfactor[i]; Kuintensity[countk]=Kintensity[i];
            countk++;
        }else{
            bool is_new=true;
            for(int j=0;j<countk;j++){
                if(fabs(Ktheta[i]-Kutheta[j])<ZERO_LIMIT){
                    is_new=false;
                    break;
                }
            }
            if(is_new){
                vector_copy(kuvector[countk], kvector[i]);
                Kutheta[countk]=Ktheta[i]; Kufactor[countk]=Kfactor[i]; Kuintensity[countk]=Kintensity[i];
                countk++;
            }

        }
    }
    for(int i=0;i<countk;i++){
        int multiplicity=0;
        for(int j=0;j<numk;j++){
            if(fabs(Kutheta[i]-Ktheta[j])<ZERO_LIMIT){
                multiplicity++;
            }
        }   
        Kuintensity[i]*=multiplicity;
    }
    for(int i=0;i<countk;i++){
        vector_copy(kvector[i], kuvector[i]);
        Ktheta[i]=Kutheta[i]; Kfactor[i]=Kufactor[i]; Kintensity[i]=Kuintensity[i];
        if(Kintensity[i]<intensity_min) intensity_min=Kintensity[i];
        if(Kintensity[i]>intensity_max) intensity_max=Kintensity[i];
    }
    numk=countk;
    deallocate_2d(kuvector, countk);
    deallocate(Kutheta); deallocate(Kufactor); deallocate(Kuintensity);
    printf("Number of unique diffraction intensity: %d\n", countk);
    printf("Range of unique diffraction intensity: %.8f %.8f", intensity_min, intensity_max);
    if(is_lorentz_flag){
        printf(" including lorentz factor\n");
    }else{
        printf("\n");
    }
}

