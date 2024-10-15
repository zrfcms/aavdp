#include "SSF.h"

SSF::SSF(RDF *rdf, double qmax, int nbin, bool is_partial)
{
    numqbin=nbin;
    numij=rdf->numij;
    qbin=qmax/(double)nbin;
	mallocate(&qij, numqbin);
	for(int i=0;i<numqbin;i++){
		qij[i]=((double)i+0.5)*qbin;
	}
    mallocate_2d(&ij, numij, numqbin);
    for(int i=0;i<numij;i++){
        for(int j=0;j<numqbin;j++){
            ij[i][j]=rdf->ij[i][j];
        }
    }
    callocate(&Smax, numij, -1.0e8);
    callocate(&Smin, numij, 1.0e8);
    callocate_2d(&Sij, numij, numqbin, 0.0);
    double constS=FOUR_PI*rdf->rhoij[0];
    for(int i=0;i<numij;i++){
        for(int j=0;j<numqbin;j++){
            for(int k=0;k<rdf->numrbin;k++){
                Sij[i][j]+=(sin(qij[j]*rdf->rij[k])/qij[j]*(rdf->gij[i][k]-1.0)*rdf->rij[k]*rdf->dr);
            }
            Sij[i][j]=1.0+constS*Sij[i][j];
            if(Smax[i]<Sij[i][j]) Smax[i]=Sij[i][j];
            if(Smin[i]>Sij[i][j]) Smin[i]=Sij[i][j];
        }
    }  
}

SSF::SSF(const char *model_path, double qmax, int nbin, bool is_partial)
{
    printf("[INFO] Starting computation of static structure factor...\n");
	QB_tools QB;
	QB_init(&QB);
	QB_read_file(&QB, model_path);
    double dq=TWO_PI/cbrt(QB.mat[0][0]*QB.mat[1][1]*QB.mat[2][2]);
    double spacingK[3]={dq, dq, dq};
    int NspacingK[3]={int(ceil(qmax/spacingK[0])), int(ceil(qmax/spacingK[1])), int(ceil(qmax/spacingK[2]))};
    printf("[INFO] Spacings along three axes in reciprocal space (Angstrom-1): %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
	
    numqbin=nbin;
	qbin=qmax/(double)nbin; 
	mallocate(&qij, numqbin);
	for(int i=0;i<numqbin;i++){
		qij[i]=((double)i+0.5)*qbin;
	}

	if(is_partial){
		numij=QB.TypeNumber*(QB.TypeNumber+1)/2+1;
	}else{
		numij=1;
	}
    callocate(&Smax, numij, -1.0e8); callocate(&Smin, numij, 1.0e8);
	callocate_2d(&ij, numij, 2, 0);
	callocate_2d(&Sij, numij, numqbin, 0.0);

	int countij=0;
	ij[countij][0]=0; ij[countij][1]=0;
	compute(&QB, qmax, qbin, spacingK, NspacingK, countij); countij++;
	if(countij<numij){
		for(int i=1;i<=QB.TypeNumber;i++){
			for(int j=i;j<=QB.TypeNumber;j++){
				ij[countij][0]=i; ij[countij][1]=j; 
				compute(&QB, qmax, qbin, spacingK, NspacingK, i, j, countij); countij++;
			}
		}
	}
    printf("[INFO] Ending computation of static structure factor\n");
}

SSF::~SSF()
{
    if(0!=numij){
        deallocate(qij);
        deallocate_2d(ij, numij);
        deallocate_2d(Sij, numij);
    }
}

void SSF::compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int pairij_id)
{
    printf("[INFO] Starting computation of static structure factor of pair *-*...\n");
    int countk=0, numk=NspacingK[0]*NspacingK[1]*NspacingK[2];
    int *nKbins; callocate(&nKbins, numqbin, 0);
    for(int i=0;i<NspacingK[0];i++){
        for(int j=0;j<NspacingK[1];j++){
            for(int k=0;k<NspacingK[2];k++){
                countk++;
                if(0==countk%10000) printf("[INFO] Completed wave vector %d of %d\n", countk, numk);
                if(0==i&&0==j&&0==k) continue;
                double K[3]={double(i)*spacingK[0], double(j)*spacingK[1], double(k)*spacingK[2]};
                double Kmag=vector_length(K);
                if(Kmag>qmax) continue;

                double cossum=0.0, sinsum=0.0;
                for(int l=0;l<QB->TotalNumber;l++){
                    double pos[3]={QB->atom[l].x, QB->atom[l].y, QB->atom[l].z};
                    double q=vector_dot(K, pos);
                    cossum+=cos(q); sinsum+=sin(q);
                }
                int ibin=int(Kmag/qbin); nKbins[ibin]++;
                Sij[pairij_id][ibin]+=(cossum*cossum+sinsum*sinsum);
            }
        }
    }
    for(int i=0;i<numqbin;i++){
        if(0!=nKbins[i]){
            Sij[pairij_id][i]/=double(nKbins[i]*QB->TotalNumber);
        }
        if(Smax[pairij_id]<Sij[pairij_id][i]) Smax[pairij_id]=Sij[pairij_id][i];
        if(Smin[pairij_id]>Sij[pairij_id][i]) Smin[pairij_id]=Sij[pairij_id][i];
    }
    deallocate(nKbins);
    printf("[INFO] Ending computation of static structure factor of pair *-*\n");
	printf("[INFO] Range of static structure factor of pair *-*: %.8f %.8f\n", Smin[pairij_id], Smax[pairij_id]);
}

void SSF::compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int typei, int typej, int pairij_id)
{
    printf("[INFO] Starting computation of static structure factor of pair %d-%d...\n", typei, typej);
    int countk=0, numk=NspacingK[0]*NspacingK[1]*NspacingK[2];
    int *nKbins; callocate(&nKbins, numqbin, 0);
    for(int i=0;i<NspacingK[0];i++){
        for(int j=0;j<NspacingK[1];j++){
            for(int k=0;k<NspacingK[2];k++){
                countk++;
                if(0==countk%10000) printf("[INFO] Completed wave vector %d of %d\n", countk, numk);
                if(0==i&&0==j&&0==k) continue;
                double K[3]={double(i)*spacingK[0], double(j)*spacingK[1], double(k)*spacingK[2]};
                double Kmag=vector_length(K);
                if(Kmag>qmax) continue;

                complex<double> sumi(0.0, 0.0), sumj(0.0, 0.0);
                for(int l=0;l<QB->TotalNumber;l++){
                    double pos[3]={QB->atom[l].x, QB->atom[l].y, QB->atom[l].z};
                    double q=vector_dot(K, pos);
                    if(QB->atom[l].type==typei) sumi+=complex<double>(cos(q), sin(q));
                    if(QB->atom[l].type==typej) sumj+=complex<double>(cos(q), -sin(q));
                }
                complex<double> sum=sumi*sumj;
                int ibin=int(Kmag/qbin); nKbins[ibin]++;
                Sij[pairij_id][ibin]+=sum.real();
            }
        }
    }
    for(int i=0;i<numqbin;i++){
        if(0!=nKbins[i]){
            Sij[pairij_id][i]/=double(nKbins[i]*QB->TotalNumber);
        }
        if(Smax[pairij_id]<Sij[pairij_id][i]) Smax[pairij_id]=Sij[pairij_id][i];
        if(Smin[pairij_id]>Sij[pairij_id][i]) Smin[pairij_id]=Sij[pairij_id][i];
    }
    deallocate(nKbins);
    printf("[INFO] Ending computation of static structure factor of pair %d-%d\n", typei, typej);
	printf("[INFO] Range of static structure factor of pair %d-%d: %.8f %.8f\n", typei, typej, Smin[pairij_id], Smax[pairij_id]);
}

// void SSF::compute(QB_tools *QB, double qmax, double qbin, double spacingK[3], int NspacingK[3], int pairij_id)
// {
//     RDF rdf(model_path, qmax, nbin, is_partial);
//     numij=rdf.numij;
//     mallocate_2d(&ij, numij, numqbin);
//     for(int i=0;i<numij;i++){
//         for(int j=0;j<numqbin;j++){
//             ij[i][j]=rdf.ij[i][j];
//         }
//     }
//     callocate_2d(&Sij, numij, numqbin, 0.0);
//     double constS=FOUR_PI*rdf.dr*rdf.rhoij[0];
//     for(int i=0;i<numij;i++){
//         for(int j=0;j<numqbin;j++){
//             for(int k=0;k<rdf.numrbin;k++){
//                 Sij[i][j]+=rdf.rij[k]*(rdf.gij[i][k]-1.0)*sin(qij[j]*rdf.rij[k])/qij[j];
//             }
//             Sij[i][j]=1.0+constS*Sij[i][j];
//         }
//     }
// }

void SSF::ssf(const char *ssf_path, const char *png_path)
{
	FILE* fp=fopen(ssf_path,"w");
	fprintf(fp,"# Static structure factor (%d data points)\n", numqbin);
	fprintf(fp,"# q\tS(q)\n");
	if(numij>1){
		fprintf(fp,"# *-*\t");
		for(int i=1;i<numij;i++){
			fprintf(fp,"%d-%d\t", ij[i][0], ij[i][1]);
		}
		fprintf(fp,"\n");
	}
	for(int i=0;i<numqbin;i++){
		fprintf(fp,"%lf\t", qij[i]);
		for(int j=0;j<numij;j++){
			fprintf(fp,"%lf\t", Sij[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
    printf("[INFO] Visualized data for static structure factor stored in %s.\n", ssf_path);
    double xlim=qbin*double(numqbin);
    double ylim=Smax[0];
    for(int i=1;i<numij;i++){
        if(ylim<Smax[i]) ylim=Smax[i];
    }
    int n_major_xtick=int(xlim/1.0), n_major_ytick=int(ylim/0.1);
    double *major_xticks, *major_yticks;
    mallocate(&major_xticks, n_major_xtick);
    mallocate(&major_yticks, n_major_ytick);
    for(int i=0;i<=n_major_xtick;i++){
        major_xticks[i]=double(i);
    }
    for(int i=0;i<=10;i++){
        major_yticks[i]=double(i)*0.1;
    }
    double height=4.0, width=6.0;
    GRAPH graph(width, height, 300);
    graph.set_xlim(0.0, xlim);
    graph.set_ylim(0.0, ylim);
    graph.set_xticks(major_xticks, n_major_xtick);
    graph.set_yticks(major_yticks, n_major_ytick);
    graph.set_tick_in(false);
    graph.set_top_visible(false);
    graph.set_right_visible(false);
    for(int i=1;i<numij;i++){
        graph.line(qij, Sij[i], numqbin);
    }
    graph.draw(png_path);
    printf("[INFO] Image for static structure factor stored in %s\n", png_path);
}

// void SSF::read(const char *nml_path)
// {
//     FILE *fp=fopen(nml_path, "r");
//     if(fp==nullptr){
//         printf("[ERROR] Unable to open file %s.\n", nml_path);
//     }
//     fseek(fp, 0, SEEK_SET);

//     char readbuff[KEY_CHAR_NUMBER]; 
//     int keyi=0, nkey=3;
//     for(int i=0;(keyi<nkey)&&(i<KEY_LINE_NUMBER);i++){
//         fscanf(fp, "%s", readbuff);
//         if('#'==readbuff[0]){
//             while(('\n'!=fgetc(fp))&&(!feof(fp)));
//             continue;
//         }
//         if(0==strcmp(readbuff, "density0")){
//             fscanf(fp, "%lf", &rho0);
//             keyi++;
//             continue;
//         }
//         if(0==strcmp(readbuff, "distance")){
//             fscanf(fp, "%lf", &dr);
//             keyi++;
//             continue;
//         }
//         if(0==strcmp(readbuff, "radial_distribution_function")){
//             fscanf(fp, "%d", &numrbin);
//             mallocate(&rij, numrbin);
//             mallocate(&gij, numrbin);
//             for(int j=0;j<numrbin;j++){
//                 fscanf(fp, "%lf", &rij[j]);
//                 fscanf(fp, "%lf", &gij[j]);
//             }
//             keyi++;
//             continue;
//         }
//     }
//     if(nkey!=keyi){
//         printf("[ERROR] Unable to recognize information from %s.\n", nml_path);
//     }
// }