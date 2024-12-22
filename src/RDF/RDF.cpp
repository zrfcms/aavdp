#include "RDF.h"

void int_to_str(char str[], int num)
{
    int i=0;
    if(num<0){
        num=-num;
        str[i++]='-';
    } 
    do{
        str[i++]=num%10+48;//0-9: 48-57
        num/=10;
    }while(num);
    str[i]='\0';
    int j=0;
    if(str[0]=='-'){
        j=1;
        i++;
    } 
    for(;j<i/2;j++){
        str[j]=str[j]+str[i-1-j];
        str[i-1-j]=str[j]-str[i-1-j];
        str[j]=str[j]-str[i-1-j];
    } 
}

RDF::RDF(char *model_path, double rmax, int nbin, bool is_partial)
{
	QB_tools QB;
	QB_init(&QB);
	QB_read_file(&QB, model_path);

	printf("[INFO] Starting computation of radial distribution function...\n");
	set_volume(QB.mat);
	numrbin=nbin;
	rbin=rmax/(double)nbin;
	callocate(&rij, numrbin, 0.0);
	for(int i=0;i<numrbin;i++){
		rij[i]=((double)i+0.5)*rbin;
	}

	if(is_partial){
		numij=QB.TypeNumber*(QB.TypeNumber+1)/2+1;
	}else{
		numij=1;
	}

	int natom=QB.TotalNumber;
	QB_pbc(&QB, rmax); QB.MCN=QB.TotalNumber;
	QB_network_init(&QB, rmax);
	callocate(&rhoij, numij, 0.0);
	callocate(&gmax, numij, 0.0);
	callocate(&gmin, numij, 1.0e8);
	callocate_2d(&ij, numij, 2, 0);
	callocate_2d(&gij, numij, numrbin, 0.0);
	int countij=0;
	compute(&QB, natom, rmax, countij);
	printf("[INFO] Density of * atoms [in Angstrom-3]: %.8f\n", rhoij[countij]);
	printf("[INFO] Range of radial distribution function of pair *-*: %.8f %.8f\n", gmin[countij], gmax[countij]);
	countij++;
	if(countij<numij){
		for(int i=1;i<=QB.TypeNumber;i++){
			for(int j=i;j<=QB.TypeNumber;j++){
				ij[countij][0]=i; ij[countij][1]=j; 
				compute(&QB, natom, rmax, i, j, countij);
				printf("[INFO] Density of j atoms of type %d [in Angstrom-3]: %.8f\n", j, rhoij[countij]);
				printf("[INFO] Range of radial distribution function of pair %d-%d: %.8f %.8f\n", i, j, gmin[countij], gmax[countij]);
				countij++;
			}
		}
	}
	QB_network_free(&QB);
	QB_pbc_clean(&QB);
	QB_free_atom(&QB);
	printf("[INFO] Ending computation of radial distribution function\n");
}

RDF::~RDF()
{
	if(0!=numij){
		deallocate(rhoij);
		deallocate(rij);
		deallocate(gmax);
		deallocate(gmin);
		deallocate_2d(ij, numij);
		deallocate_2d(gij, numij);
	}
}

void RDF::set_volume(double mat[3][3])
{
	double odd=0.0, even=0.0;
	for(int i=0;i<3;i++){
        odd+=mat[0][i]*mat[1][(i+1)%3]*mat[2][(i+2)%3];
    }
	for(int i=0;i<3;i++){
        even+=mat[0][i]*mat[1][(i+2)%3]*mat[2][(i+1)%3];
    }
	vol=fabs(odd-even);
}

void RDF::compute(QB_tools *QB, int natom, double rmax, int pairij_id)
{
    printf("[INFO] Starting computation of radial distribution function of pair *-*...\n");
	clock_t start, finish;
    start=clock();
	int *natombin; callocate(&natombin, numrbin, 0);
	for(int i=0;i<natom;i++){
		QB_list_build(QB, i, rmax);
		double posi[3]={QB->atom[i].x, QB->atom[i].y, QB->atom[i].z};
		for(int j=0;j<QB->neb.num;j++){
			int    id=QB->neb.id[j];
			double posj[3]={QB->atom[id].x, QB->atom[id].y, QB->atom[id].z};
			double dist[3];
			vector_difference(dist, posi, posj);
			id=int(vector_length(dist)/rbin);
			if(id>=0&&id<numrbin){
				natombin[id]++;
			}
		}
		QB_list_clear(QB);		
	}
	finish=clock();
    printf("[INFO] Ending computation of radial distribution function of pair *-*\n");
	printf("[INFO] Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);

	rhoij[pairij_id]=(double)natom/vol;
	double constg=THREE_QUARTER_PIINV/pow(rbin, 3)/rhoij[pairij_id]/(double)natom;
	for(int i=0;i<numrbin;i++){
		gij[0][i]=constg*natombin[i]/double(pow(i+1, 3)-pow(i, 3));
		if(gmax[pairij_id]<gij[0][i]) gmax[pairij_id]=gij[pairij_id][i];
		if(gmin[pairij_id]>gij[0][i]) gmin[pairij_id]=gij[pairij_id][i];
	}
	deallocate(natombin);
}

void RDF::compute(QB_tools *QB, int natom, double rmax, int typei, int typej, int pairij_id)
{
    printf("[INFO] Starting computation of radial distribution function of pair %d-%d...\n", typei, typej);
	clock_t start, finish;
    start=clock();
	int natomi=0, natomj=0;
	int *natombin; callocate(&natombin, numrbin, 0);
	for(int i=0;i<natom;i++){
		if(QB->atom[i].type==typei){
			natomi++;
			QB_list_build(QB, i, rmax);
			double posi[3]={QB->atom[i].x, QB->atom[i].y, QB->atom[i].z};
			for(int j=0;j<QB->neb.num;j++){
				int    id=QB->neb.id[j];
				if(QB->atom[id].type==typej){
					double posj[3]={QB->atom[id].x, QB->atom[id].y, QB->atom[id].z};
					double dist[3];
					vector_difference(dist, posi, posj);
					id=vector_length(dist)/rbin;
					if(id>=0&&id<numrbin){
						natombin[id]++;
					}
				}
			}
			QB_list_clear(QB);
		}else if(QB->atom[i].type==typej){
			natomj++;
		}	
	}
	if(typei==typej) natomj=natomi;
	finish=clock();
    printf("[INFO] Ending computation of radial distribution function of pair %d-%d\n", typei, typej);
	printf("[INFO] Computation time [s]: %.8f\n", double(finish-start)/CLOCKS_PER_SEC);

	rhoij[pairij_id]=natomj/vol;
	double constg=THREE_QUARTER_PIINV/pow(rbin, 3)/rhoij[pairij_id]/(double)natomi;
	for(int i=0;i<numrbin;i++){
		gij[pairij_id][i]=constg*natombin[i]/double(pow(i+1, 3)-pow(i, 3));
		if(gmax[pairij_id]<gij[pairij_id][i]) gmax[pairij_id]=gij[pairij_id][i];
		if(gmin[pairij_id]>gij[pairij_id][i]) gmin[pairij_id]=gij[pairij_id][i];
	}
	deallocate(natombin);
}

void RDF::img(char *png_path, double *x, double *y, int num, double ymin, double ymax)
{
    double xmin=0.0, xmax=ceil(x[num-1]);
    ymin=floor(ymin); ymax=ceil(ymax);
    int n_major_xtick=int(xmax-xmin)+1, n_major_ytick=int(ymax-ymin)+1;
    double *major_xticks, *major_yticks;
    callocate(&major_xticks, n_major_xtick, 0.0);
    callocate(&major_yticks, n_major_ytick, 0.0);
    for(int i=0;i<n_major_xtick;i++){
        major_xticks[i]=xmin+double(i);
    }
    for(int i=0;i<n_major_ytick;i++){
        major_yticks[i]=ymin+double(i);
    }
    double height=4.0, width=6.0;
    GRAPH graph(width, height, 300);
    graph.set_xlim(major_xticks[0], major_xticks[n_major_xtick-1]);
    graph.set_ylim(major_yticks[0], major_yticks[n_major_ytick-1]);
    graph.set_xticks(major_xticks, n_major_xtick);
    graph.set_yticks(major_yticks, n_major_ytick);
    graph.set_tick_in(false);
    graph.set_top_visible(false);
    graph.set_right_visible(false);
    graph.line(x, y, num);
    graph.draw(png_path);
	deallocate(major_xticks);
	deallocate(major_yticks);
}

void RDF::rdf(char *rdf_path)
{
	FILE* fp=fopen(rdf_path,"w");
	fprintf(fp,"# r\tg(r) (%d points)\n", numrbin);
	if(numij>1){
		fprintf(fp,"# *-*\t");
		for(int i=1;i<numij;i++){
			fprintf(fp,"%d-%d\t", ij[i][0], ij[i][1]);
		}
		fprintf(fp,"\n");
	}
	for(int i=0;i<numrbin;i++){
		fprintf(fp,"%.8f\t", rij[i]);
		for(int j=0;j<numij;j++){
			fprintf(fp,"%.8f\t", gij[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
    printf("[INFO] Information for radial distribution function stored in %s\n", rdf_path);

    char png_path[strlen(rdf_path)+20];
	double ggmin=1.0e8, ggmax=0.0;
	for(int i=0;i<numij;i++){
		if(ggmin>gmin[i]) ggmin=gmin[i];
		if(ggmax<gmax[i]) ggmax=gmax[i];
	}
	for(int i=0;i<numij;i++){
		int typei=ij[i][0], typej=ij[i][1];
		strcpy(png_path, rdf_path);
		char stri[5], strj[5]; int_to_str(stri, typei); int_to_str(strj, typej); 
		strcat(png_path, "."); strcat(png_path, stri); strcat(png_path, "-"); strcat(png_path, strj);
		strcat(png_path, ".png");
		img(png_path, rij, gij[i], numrbin, ggmin, ggmax);
		printf("[INFO] Image for radial distribution function stored in %s\n", png_path);
	}
}

SSF::SSF(RDF *rdf, double qmax, int nbin, bool is_partial)
{
	printf("[INFO] Starting computation of static structure factor...\n");
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
                Sij[i][j]+=(sin(qij[j]*rdf->rij[k])/qij[j]*(rdf->gij[i][k]-1.0)*rdf->rij[k]*rdf->rbin);
            }
            Sij[i][j]=1.0+constS*Sij[i][j];
            if(Smax[i]<Sij[i][j]) Smax[i]=Sij[i][j];
            if(Smin[i]>Sij[i][j]) Smin[i]=Sij[i][j];
        }
    }
	printf("[INFO] Range of static structure factor of pair *-*: %.8f %.8f\n", Smin[0], Smax[0]);
	for(int i=1;i<numij;i++){
		printf("[INFO] Range of static structure factor of pair %d-%d: %.8f %.8f\n", ij[i][0], ij[i][1], Smin[i], Smax[i]);
	}
	printf("[INFO] Ending computation of static structure factor\n");  
}

SSF::SSF(char *model_path, double qmax, int nbin, bool is_partial)
{
	QB_tools QB;
	QB_init(&QB);
	QB_read_file(&QB, model_path);
	printf("[INFO] Starting computation of static structure factor...\n");
    double dq=TWO_PI/cbrt(QB.mat[0][0]*QB.mat[1][1]*QB.mat[2][2]);
    double spacingK[3]={dq, dq, dq};
    int NspacingK[3]={int(ceil(qmax/spacingK[0])), int(ceil(qmax/spacingK[1])), int(ceil(qmax/spacingK[2]))};
    printf("[INFO] Spacings along a*, b*, and c* in reciprocal space [Angstrom-1]: %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
	
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
	compute(&QB, qmax, spacingK, NspacingK, countij);
	printf("[INFO] Range of static structure factor of pair *-*: %.8f %.8f\n", Smin[countij], Smax[countij]);
	countij++;
	if(countij<numij){
		for(int i=1;i<=QB.TypeNumber;i++){
			for(int j=i;j<=QB.TypeNumber;j++){
				ij[countij][0]=i; ij[countij][1]=j; 
				compute(&QB, qmax, spacingK, NspacingK, i, j, countij);
				printf("[INFO] Range of static structure factor of pair %d-%d: %.8f %.8f\n", i, j, Smin[countij], Smax[countij]);
				countij++;
			}
		}
	}
    printf("[INFO] Ending computation of static structure factor\n");
}

SSF::~SSF()
{
    if(0!=numij){
        deallocate(qij);
		deallocate(Smax);
		deallocate(Smin);
        deallocate_2d(ij, numij);
        deallocate_2d(Sij, numij);
    }
}

void SSF::compute(QB_tools *QB, double qmax, double spacingK[3], int NspacingK[3], int pairij_id)
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
}

void SSF::compute(QB_tools *QB, double qmax, double spacingK[3], int NspacingK[3], int typei, int typej, int pairij_id)
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
}

void SSF::img(char *png_path, double *x, double *y, int num, double ymin, double ymax)
{
    double xmin=0.0, xmax=ceil(x[num-1]);
    ymin=floor(ymin); ymax=ceil(ymax);
    int n_major_xtick=int(xmax-xmin)+1, n_major_ytick=int(ymax-ymin)+1;
    double *major_xticks, *major_yticks;
    callocate(&major_xticks, n_major_xtick, 0.0);
    callocate(&major_yticks, n_major_ytick, 0.0);
    for(int i=0;i<n_major_xtick;i++){
        major_xticks[i]=xmin+double(i);
    }
    for(int i=0;i<n_major_ytick;i++){
        major_yticks[i]=ymin+double(i);
    }
    double height=4.0, width=6.0;
    GRAPH graph(width, height, 300);
    graph.set_xlim(major_xticks[0], major_xticks[n_major_xtick-1]);
    graph.set_ylim(major_yticks[0], major_yticks[n_major_ytick-1]);
    graph.set_xticks(major_xticks, n_major_xtick);
    graph.set_yticks(major_yticks, n_major_ytick);
    graph.set_tick_in(false);
    graph.set_top_visible(false);
    graph.set_right_visible(false);
    graph.line(x, y, num);
    graph.draw(png_path);
	deallocate(major_xticks);
	deallocate(major_yticks);
}

void SSF::ssf(char *ssf_path)
{
	FILE* fp=fopen(ssf_path,"w");
	fprintf(fp,"# q\tS(q) (%d points)\n", numqbin);
	if(numij>1){
		fprintf(fp,"# *-*\t");
		for(int i=1;i<numij;i++){
			fprintf(fp,"%d-%d\t", ij[i][0], ij[i][1]);
		}
		fprintf(fp,"\n");
	}
	for(int i=0;i<numqbin;i++){
		fprintf(fp,"%.8f\t", qij[i]);
		for(int j=0;j<numij;j++){
			fprintf(fp,"%.8f\t", Sij[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
    printf("[INFO] Information for static structure factor stored in %s\n", ssf_path);

    char png_path[strlen(ssf_path)+20];
	double SSmin=1.0e8, SSmax=-1.0e8;
	for(int i=0;i<numij;i++){
		if(SSmin>Smin[i]) SSmin=Smin[i];
		if(SSmax<Smax[i]) SSmax=Smax[i];
	}
	for(int i=0;i<numij;i++){
		int typei=ij[i][0], typej=ij[i][1];
		strcpy(png_path, ssf_path);
		char stri[5], strj[5]; int_to_str(stri, typei); int_to_str(strj, typej); 
		strcat(png_path, "."); strcat(png_path, stri); strcat(png_path, "-"); strcat(png_path, strj);
		strcat(png_path, ".png");
		img(png_path, qij, Sij[i], numqbin, SSmin, SSmax);
		printf("[INFO] Image for static structure factor stored in %s\n", png_path);
	}
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