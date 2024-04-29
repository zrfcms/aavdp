#include "SSF.h"

SSF::SSF(const char *nml_path, double qmax, double nbin)
{
    read(nml_path);
    double dq=qmax/double(nbin);
    numqbin=nbin;
    callocate(&qij, numqbin, 0.0);
    callocate(&Sij, numqbin, 0.0);
    double constS=4*PI*rho0*dr;
    for(int i=0;i<numqbin;i++){
        qij[i]=dq*((double)i+0.5);
        for(int j=0;j<numrbin;j++){
            Sij[i]+=rij[j]*rij[j]*(gij[j]-1.0)*sin(qij[i]*rij[j])/(qij[i]*rij[j]);
        }
        Sij[i]*=constS;
    }
}

SSF::~SSF()
{
    if(0!=numrbin){
        deallocate(rij);
        deallocate(gij);
    }
    if(0!=numqbin){
        deallocate(qij);
        deallocate(Sij);
    }
}

void SSF::read(const char *nml_path)
{
    FILE *fp=fopen(nml_path, "r");
    if(fp==nullptr){
        printf("[ERROR] Unable to open file %s.\n", nml_path);
    }
    fseek(fp, 0, SEEK_SET);

    char readbuff[KEY_CHAR_NUMBER]; 
    int keyi=0, nkey=3;
    for(int i=0;(keyi<nkey)&&(i<KEY_LINE_NUMBER);i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, "density0")){
            fscanf(fp, "%lf", &rho0);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, "distance")){
            fscanf(fp, "%lf", &dr);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, "radial_distribution_function")){
            fscanf(fp, "%d", &numrbin);
            mallocate(&rij, numrbin);
            mallocate(&gij, numrbin);
            for(int j=0;j<numrbin;j++){
                fscanf(fp, "%lf", &rij[j]);
                fscanf(fp, "%lf", &gij[j]);
            }
            keyi++;
            continue;
        }
    }
    if(nkey!=keyi){
        printf("[ERROR] Unable to recognize information from %s.\n", nml_path);
    }
}

void SSF::ssf(const char *ssf_path)
{
	FILE* f=fopen(ssf_path,"w");
	fprintf(f,"# Static Structure Factor (%d data points)\n", numqbin);
	fprintf(f,"# q\tS(q)\n");
	for(int i=0;i<numqbin;i++){
		fprintf(f,"%lf\t%lf\t", qij[i], Sij[i]);
		fprintf(f,"\n");
	}
	fclose(f);
}