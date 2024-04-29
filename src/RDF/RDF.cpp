#include "RDF.h"

RDF::RDF(const char *model_path, double rmax, int nbin, bool is_partial_flag)
{
	QB_tools QB;
	QB_init(&QB);
	QB_read_file(&QB, model_path);

    printf("Starting...\n");
	dr=rmax/(double)nbin; numrbin=nbin;
	volume(&QB);
	mallocate(&rij, numrbin);
	for(int i=0;i<numrbin;i++){
		rij[i]=((double)i+0.5)*dr;
	}
	if(is_partial_flag){
		numij=QB.TypeNumber*(QB.TypeNumber+1)/2;
		mallocate_2d(&ij, numij, 2);
		mallocate(&rho0ij, numij);
		mallocate_2d(&gij, numij, numrbin);
		int countij=0;
		for(int i=1;i<=QB.TypeNumber;i++){
			for(int j=i;j<=QB.TypeNumber;j++){
				ij[countij][0]=i; ij[countij][1]=j;
				compute(&QB, rmax, i, j, countij);
				countij++;
			}
		}
	}else{
		numij=1;
		callocate_2d(&ij, 1, 2, 0);
		callocate(&rho0ij, 1, 0.0);
		mallocate_2d(&gij, 1, numrbin);
		compute(&QB, rmax);
	}
	printf("Ending...\n");
	QB_free_atom(&QB);
}

RDF::~RDF()
{
	if(0!=numij){
		deallocate_2d(ij, numij);
		deallocate(rho0ij);
		deallocate(rij);
		deallocate_2d(gij, numij);
	}
}

void RDF::volume(QB_tools *QB)
{
	double odd=0.0, even=0.0;
	for(int i=0;i<3;i++){
        odd+=QB->mat[0][i]*QB->mat[1][(i+1)%3]*QB->mat[2][(i+2)%3];
    }
	for(int i=0;i<3;i++){
        even+=QB->mat[0][i]*QB->mat[1][(i+2)%3]*QB->mat[2][(i+1)%3];
    }
	vol=fabs(odd-even);
}

void RDF::compute(QB_tools *QB, double rmax, int typei, int typej, int pair_id)
{
	int natom=QB->TotalNumber;
	QB_pbc(QB, rmax);
	QB->MCN=QB->TotalNumber;
	QB_network_init(QB, rmax);
	int natomi=0, natomj=0;
	int *natombin; callocate(&natombin, numrbin, 0);
    printf("Starting computation of pair %d-%d...\n", typei, typej);
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
					id=vector_length(dist)/dr;
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
	QB_network_free(QB);
	QB_pbc_clean(QB);

	rho0ij[pair_id]=natomi/vol;
	double constg=3.0/(4.0*PI)/pow(dr, 3)/rho0ij[pair_id]/natomj;
	double gmax=0.0, gmin=1.0e8;
	for(int i=0;i<numrbin;i++){
		gij[pair_id][i]=constg*natombin[i]/(pow(double(i+1),3)-pow(double(i),3));
		if(gmax<gij[pair_id][i]) gmax=gij[pair_id][i];
		if(gmin>gij[pair_id][i]) gmin=gij[pair_id][i];
	}
	deallocate(natombin);
    printf("Ending computation of pair %d-%d\n", typei, typej);
	printf("Density of centering atoms %d: %.8f\n", typei, rho0ij[pair_id]);
	printf("Range of radial distribution function of pair %d-%d: %.8f %.8f\n", typei, typej, gmin, gmax);
}

void RDF::compute(QB_tools *QB, double rmax)
{
	int natom=QB->TotalNumber;
	QB_pbc(QB, rmax);
	QB->MCN=QB->TotalNumber;
	QB_network_init(QB, rmax);
	int natomi=0, natomj=0;
	int *natombin; callocate(&natombin, numrbin, 0);
    printf("Starting computation of pair *-*...\n");
	for(int i=0;i<natom;i++){
		QB_list_build(QB, i, rmax);
		double posi[3]={QB->atom[i].x, QB->atom[i].y, QB->atom[i].z};
		for(int j=0;j<QB->neb.num;j++){
			int    id=QB->neb.id[j];
			double posj[3]={QB->atom[id].x, QB->atom[id].y, QB->atom[id].z};
			double dist[3];
			vector_difference(dist, posi, posj);
			id=int(vector_length(dist)/dr);
			if(id>=0&&id<numrbin){
				natombin[id]++;
			}
		}
		QB_list_clear(QB);		
	}
	QB_network_free(QB);
	QB_pbc_clean(QB);

	rho0ij[0]=(double)natom/vol;
	double constg=3.0/(4.0*PI)/pow(dr, 3)/rho0ij[0]/(double)natom;
	double gmax=0.0, gmin=1.0e8;
	for(int i=0;i<numrbin;i++){
		gij[0][i]=constg*natombin[i]/(pow(double(i+1),3)-pow(double(i),3));
		if(gmax<gij[0][i]) gmax=gij[0][i];
		if(gmin>gij[0][i]) gmin=gij[0][i];
	}
	deallocate(natombin);
    printf("Ending computation of pair *-*\n");
	printf("Density of centering atoms: %.8f\n", rho0ij[0]);
	printf("Range of radial distribution function of pair *-*: %.8f %.8f\n", gmin, gmax);
}

void RDF::rdf(const char *rdf_path)
{
	FILE* f=fopen(rdf_path,"w");
	fprintf(f,"# Radial Distribution Function (%d data points)\n", numrbin);
	fprintf(f,"# r\tg(r)\n");
	if(numij>1){
		fprintf(f,"#");
		for(int i=0;i<numij;i++){
			fprintf(f," %d-%d\t", ij[i][0], ij[i][1]);
		}
		fprintf(f,"\n");
	}
	for(int i=0;i<numrbin;i++){
		fprintf(f,"%lf\t", rij[i]);
		for(int j=0;j<numij;j++){
			fprintf(f,"%lf\t", gij[j][i]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
    printf("Visualized data stored in %s.\n", rdf_path);
}

void RDF::nml(const char *nml_path)
{
	FILE* f=fopen(nml_path, "w");
	fprintf(f,"distance %lf\n", dr);
	fprintf(f,"density0 %lf\n", rho0ij[0]);
	fprintf(f,"radial_distribution_function %d\n", numrbin);
	for(int i=0;i<numrbin;i++){
		fprintf(f,"%lf\t%lf\n", rij[i], gij[0][i]);
	}
	fclose(f);
    printf("Information for static structure factor stored in %s.\n", nml_path);
}