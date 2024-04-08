#include "EBSD_MATH.h"

bool is_in_vector(const int value, const int vec[], const int num)
{
    for(int i=0;i<num;i++){
        if(value==vec[i]) return true;
    }
    return false;
}

template <typename T>
void normalize_vector(T n_v[], T v[], int num)
{
    T sum_v2, len_v;
    sum_v2=0.0;
    for(int i=0;i<num;i++){
        sum_v2+=v[i]*v[i];
    }
    len_v=sqrt(sum_v2);
    if(0.0!=len_v){
        n_v[0]=v[0]/len_v; n_v[1]=v[1]/len_v; n_v[2]=v[2]/len_v;
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}

void compute_square_Lambert(double xy[2], int &ierr, double xyz[3]) 
{
    ierr=0;
    xy[0]=0.0; xy[1]=0.0;
    if(fabs(1.0-(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]))>1e-12){ //Check whether the xyz lies on the unit sphere
        ierr=1;
    }else{
        if(1.0!=fabs(xyz[2])){
            double q;
            if(fabs(xyz[1])<=fabs(xyz[0])){
                q=fabs(xyz[0])/xyz[0]*sqrt(2.0*(1.0-fabs(xyz[2])));
                xy[0]=q*PI_SQRT_HALF; 
                xy[1]=q*atan(xyz[1]/xyz[0])/PI_SQRT_HALF;
            }else{
                q=fabs(xyz[1])/xyz[1]*sqrt(2.0*(1.0-fabs(xyz[2])));
                xy[0]=q*atan(xyz[0]/xyz[1])/PI_SQRT_HALF; 
                xy[1]=q*PI_SQRT_HALF;
            }
        }
    }
    xy[0]=xy[0]/PI_HALF_SQRT; xy[1]=xy[1]/PI_HALF_SQRT;
}

void compute_hexagonal_Lambert(double xy[2], int &ierr, double xyz[3])
{
    ierr=0;
    double cxy[2]={0.0, 0.0};
    if(fabs(1.0-(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]))>1e-12){
        ierr=1;
    }else{
        if(1.0!=fabs(xyz[2])){
            double q0=sqrt(2.0/(1.0+fabs(xyz[2]))), q;
            double x0=q0*xyz[1]+1.0e-4, y0=q0*xyz[0]+1.0e-4;
            double len=sqrt(x0*x0+y0*y0);
            if(x0<0.0){
                len=-len;
            }
            int sextant=get_sextant(x0, y0);
            double x, y;
            switch(sextant)
            {
            case 0:
            case 3:
                q0=PI_PREE*len;
                x=PI_HALF_SQRT*q0;
                if(0.0==x0){
                    y=PI_PREF*PI*0.5;
                }else{
                    y=PI_PREF*q0*atan(y0/x0);
                }
                break;
            case 1:
            case 4:
                q0=PI_PREA*len;
                q=atan((y0-SQRT_3*x0)/(x0+SQRT_3*y0));
                x=SQRT_3*q0*(PI/6.0-q); 
                y=q0*(0.5*PI+q);
                break;
            case 2:
            case 5:
                q0=PI_PREA*len;
                q=atan((y0+SQRT_3*x0)/(x0-SQRT_3*y0));
                x=SQRT_3*q0*(PI/6.0+q); 
                y=q0*(-0.5*PI+q);
                break;
            default:
                printf("Error! Unrecognized sextant %d (not between 0 and 5).", sextant);
                exit(EXIT_FAILURE);
            }
            cxy[0]=y; 
            cxy[1]=x;
        }
    }
    xy[0]=(cxy[0]+cxy[1]*SQRT_3_INVERSE)/PI_PREG;
    xy[1]=cxy[1]*2.0*SQRT_3_INVERSE/PI_PREG;
}

void compute_Lambert_interpolation(double xyz[3], int nump, bool hexagonal_flag, int &ix, int &iy, int &ixp, int &iyp, double &dx, double &dy, double &dxm, double &dym)
{
    double cxyz[3]={xyz[0], xyz[1], xyz[2]};
    if(cxyz[2]<0.0){
        cxyz[0]=-cxyz[0]; cxyz[1]=-cxyz[1]; cxyz[2]=-cxyz[2];
    }
    normalize_vector(cxyz, cxyz);
    double xy[2]; int ierr;
    if(hexagonal_flag){
        compute_hexagonal_Lambert(xy, ierr, cxyz);
    }else{
        compute_square_Lambert(xy, ierr, cxyz);
    }
    if(0!=ierr){
        printf("Error! Unable to compute square Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
        exit(EXIT_FAILURE);
    }
    xy[0]*=nump; xy[1]*=nump; 
    ix=int(nump+xy[0])-nump; iy=int(nump+xy[1])-nump;
    ixp=ix+1; iyp=iy+1;
    if(ixp>nump) ixp=ix;
    if(iyp>nump) iyp=iy;
    if(ix<-nump) ix=ixp;
    if(iy<-nump) iy=iyp;
    dx=xy[0]-ix; dy=xy[1]-iy; 
    dxm=1.0-dx; dym=1.0-dy;
    ixp+=nump; iyp+=nump; ix+=nump; iy+=nump;
}

int get_sextant(double x, double y)
{
    double xx=fabs(SQRT_3_INVERSE*x);
    if(x>=0.0){
        if(fabs(y)<=xx){
            return 0;
        }else{
            if(y>xx){
                return 1;
            }else{
                return 5;
            }
        }
    }else{
        if(fabs(y)<=xx){
            return 3;
        }else{
            if(y>xx){
                return 2;
            }else{
                return 4;
            }
        }
    }
    return -1;
}

void compute_sphere_from_square_Lambert(double xyz[3], int &ierr, double xy[2])
{
    ierr=0;
    xyz[0]=xyz[1]=xyz[2]=0.0;
    double cxy[2]={xy[0]*PI_HALF_SQRT, xy[1]*PI_HALF_SQRT};
    //make sure that the input point lies inside the square with certain length
    if(fmax(fabs(cxy[0]), fabs(cxy[1]))>PI_HALF_SQRT){
        ierr=1;
    }else{
        if(0.0==fmax(fabs(cxy[0]), fabs(cxy[1]))){
            xyz[2]=1.0;
        }else{
            double x, y;
            if(fabs(cxy[0])<=fabs(cxy[1])){
                x=2.0*PI_INVERSE*cxy[1]*sqrt(PI-cxy[1]*cxy[1]);
                y=0.25*PI*xy[0]/xy[1];
                xyz[0]=x*sin(y); xyz[1]=x*cos(y); xyz[2]=1.0-2.0*PI_INVERSE*cxy[1]*cxy[1];
            }else{
                x=2.0*PI_INVERSE*cxy[0]*sqrt(PI-cxy[0]*cxy[0]);
                y=0.25*PI*xy[1]/xy[0];
                xyz[0]=x*cos(y); xyz[1]=x*sin(y); xyz[2]=1.0-2.0*PI_INVERSE*cxy[0]*cxy[0];
            }
            normalize_vector(xyz, xyz);
        }
    }
}

void compute_sphere_from_hexagonal_Lambert(double xyz[3], int &ierr, double xy[2])
{
    ierr=0;
    double cxy[2]={PI_PREG*(xy[0]-0.5*xy[1]), SQRT_HALF_3*PI_PREG*xy[1]};
    if(0.0==fmax(fabs(cxy[0]), fabs(cxy[1]))){
        xyz[0]=xyz[1]=0.0; 
        xyz[2]=1.0;
    }else{
        double cyx[2]={cxy[1], cxy[0]};
        int sextant=get_sextant(cyx[0], cyx[1]);
        double q1, q2;
        double x, y;
        switch(sextant)
        {
        case 0:
        case 3:
            q1=PI_PREC*cyx[1]/cyx[0];
            x=PI_PREB*cyx[0]*cos(q1); 
            y=PI_PREB*cyx[0]*sin(q1);
            break;
        case 1:
        case 4:
            q1=cyx[0]+SQRT_3*cyx[1];
            q2=PI_PRED*cyx[0]/q1;
            x=PI_PREA*q1*sin(q2); 
            y=PI_PREA*q1*cos(q2);
            break;
        case 2:
        case 5:
            q1=cyx[0]-SQRT_3*cyx[1];
            q2=PI_PRED*cyx[0]/q1;
            x=PI_PREA*q1*sin(q2); 
            y=-PI_PREA*q1*cos(q2);
            break;
        default:
            printf("Error! Unrecognized sextant %d (not between 0 and 5).", sextant);
            exit(EXIT_FAILURE);
        }
        double q=x*x+y*y;
        //make sure that the input point lies inside the hexagon
        if(q>4.0){
            xyz[0]=xyz[1]=xyz[2]=0.0;
            ierr=1;
        }else{
            xyz[0]=0.5*x*sqrt(4.0-q); 
            xyz[1]=0.5*y*sqrt(4.0-q); 
            xyz[2]=1.0-0.5*q; 
        }
        double temp=xyz[0]; 
        xyz[0]=xyz[1]; 
        xyz[1]=temp;
    }
}

void compute_sphere_from_stereographic_projection(double xyz[3], int &ierr, double xy[2], double radius)
{
    ierr=0;
    xyz[0]=0.0; xyz[1]=0.0; xyz[2]=0.0;
    if(0.0==fmax(fabs(xy[0]), fabs(xy[1]))){
        xyz[2]=1.0; 
    }else{
        double sum_q2=xy[0]*xy[0]+xy[1]*xy[1];
        if(sum_q2>radius*radius){
            ierr=1;
        }else{
            xyz[0]=2.0*xy[0]; xyz[1]=2.0*xy[1]; xyz[2]=1.0-sum_q2;
            for(int i=0;i<3;i++){
                xyz[i]=xyz[i]/(radius*radius+sum_q2);
            }
        }
    }
}

void create_image(const char* png_path, unsigned char *pixels, int width, int height)
{	
    png_structp png_ptr;  
    png_infop info_ptr;  
    FILE *fp=fopen(png_path, "wb");  
    if(!fp){
        printf("Error! Unable to create %s using fopen.\n", png_path);
        exit(EXIT_FAILURE);
    }
    png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);  
    if(png_ptr==NULL){  
        printf("Error! Unable to create %s using png_create_write_struct.\n", png_path);
        fclose(fp);
        exit(EXIT_FAILURE);
    }  
    info_ptr=png_create_info_struct(png_ptr);
    if(info_ptr==NULL){  
        printf("Error! Unable to create %s using png_create_info_struct.\n", png_path);
        png_destroy_write_struct(&png_ptr, NULL);  
        exit(EXIT_FAILURE);
    }  
    png_init_io(png_ptr, fp);  
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, 
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE); 
    png_colorp palette=(png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH*sizeof(png_color));
    if(!palette){
        printf("Error! Unable to create palette using png_malloc.\n");
        fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        exit(EXIT_FAILURE);
    }
    png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);  
    png_write_info(png_ptr, info_ptr);  
    png_set_packing(png_ptr);
    png_bytepp rows=(png_bytepp)png_malloc(png_ptr, height*sizeof(png_bytep));
    for(int i=0;i<height;++i){
        rows[height-i-1]=(png_bytep)(pixels+(i)*width*3);
    }
    png_write_image(png_ptr, rows);  
    delete[] rows;
    png_write_end(png_ptr, info_ptr);  
    png_free(png_ptr, palette);  
    palette=NULL;  
    png_destroy_write_struct(&png_ptr, &info_ptr);  
    fclose(fp);   
}

