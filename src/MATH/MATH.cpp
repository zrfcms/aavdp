#include "MATH.h"

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

void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num)
{
    for(int i=0;i<num;i++){
        strcat(file_path, exts[i]);
    }
}

void split_path(char name[], char ext[], const char file_path[])
{
    char *ch=strrchr(file_path,'.');
    strcpy(name, file_path);
    strcpy(ext, ch);
    if(ch!=nullptr){
        int pos=ch-file_path;
        name[pos]='\0';
    }
}

void compute_square_Lambert(double xy[2], int &ierr, double xyz[3]) 
{
    ierr=0;
    xy[0]=0.0; xy[1]=0.0;
    if(fabs(1.0-vector_dot(xyz, xyz))>1.0e-8){ //Check whether the xyz lies on the unit sphere
        ierr=1;
    }else{
        if(fabs(1.0-fabs(xyz[2]))>1.0e-8){
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

void compute_sphere_from_square_Lambert(double xyz[3], int &ierr, double xy[2])
{
    ierr=0;
    xyz[0]=xyz[1]=xyz[2]=0.0;
    double cxy[2]={xy[0]*PI_HALF_SQRT, xy[1]*PI_HALF_SQRT};
    if(fmax(fabs(cxy[0]), fabs(cxy[1]))>PI_HALF_SQRT){//make sure that the input point lies inside the square with half length PI_HALF_SQRT
        ierr=1;
    }else{
        if(fmax(fabs(cxy[0]), fabs(cxy[1]))<1.0e-8){
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
            vector_normalize(xyz, xyz);
        }
    }
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

void compute_hexagonal_Lambert(double xy[2], int &ierr, double xyz[3])
{
    ierr=0;
    double x=0.0, y=0.0;
    if(fabs(1.0-vector_dot(xyz, xyz))>1.0e-8){
        ierr=1;
    }else{
        if(fabs(fabs(xyz[2])-1.0)>1.0e-8){
            double q0=sqrt(2.0/(1.0+fabs(xyz[2])));
            double x0=q0*xyz[1]+1.0e-4, y0=q0*xyz[0]+1.0e-4;
            int sextant=get_sextant(x0, y0);
            double q=sqrt(x0*x0+y0*y0);
            if(x0<0.0) q=-q;
            switch(sextant)
            {
            case 0:
            case 3:
                q=PI_PREE*q;
                x=PI_HALF_SQRT*q;
                if(fabs(x0)<1.0e-8){
                    y=PI_PREF*PI*0.5;
                }else{
                    y=PI_PREF*q*atan(y0/x0);
                }
                break;
            case 1:
            case 4:
                q=PI_PREA*q;
                q0=atan((y0-SQRT_3*x0)/(x0+SQRT_3*y0));
                x=SQRT_3*q*(PI/6.0-q0); 
                y=q*(0.5*PI+q0);
                break;
            case 2:
            case 5:
                q=PI_PREA*q;
                q0=atan((y0+SQRT_3*x0)/(x0-SQRT_3*y0));
                x=SQRT_3*q*(PI/6.0+q0); 
                y=q*(-0.5*PI+q0);
                break;
            default:
                printf("[ERROR] Unrecognized hexagonal sextant %d (not between 0 and 5).", sextant);
                exit(EXIT_FAILURE);
            }
        }
    }
    xy[0]=(y+x*SQRT_3_INVERSE)/PI_PREG;
    xy[1]=x*2.0*SQRT_3_INVERSE/PI_PREG;
}

void compute_sphere_from_hexagonal_Lambert(double xyz[3], int &ierr, double xy[2])
{
    ierr=0;
    //double x0=PI_PREG*SQRT_HALF_3*xy[0], y0=PI_PREG*(xy[1]-0.5*xy[0]);
    double x0=PI_PREG*SQRT_HALF_3*xy[1], y0=PI_PREG*(xy[0]-0.5*xy[1]);
    if(fmax(fabs(x0), fabs(y0))<1.0e-8){
        xyz[0]=xyz[1]=0.0; xyz[2]=1.0;
    }else{
        int sextant=get_sextant(x0, y0);
        double q, q0;
        double x, y;
        switch(sextant)
        {
        case 0:
        case 3:
            q=x0;
            q0=PI_PREC*y0/q;
            x=PI_PREB*q*cos(q0);
            y=PI_PREB*q*sin(q0);
            break;
        case 1:
        case 4:
            q=x0+SQRT_3*y0;
            q0=PI_PRED*x0/q;
            x=PI_PREA*q*sin(q0); 
            y=PI_PREA*q*cos(q0);
            break;
        case 2:
        case 5:
            q=x0-SQRT_3*y0;
            q0=PI_PRED*x0/q;
            x=PI_PREA*q*sin(q0); 
            y=-PI_PREA*q*cos(q0);
            break;
        default:
            printf("[ERROR] Unrecognized hexagonal sextant %d (not between 0 and 5).", sextant);
            exit(EXIT_FAILURE);
        }
        q=x*x+y*y;
        if(q>4.0){//make sure that the input point lies inside the hexagon
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
    if(fmax(fabs(xy[0]), fabs(xy[1]))<1.0e-8){
        xyz[2]=1.0; 
    }else{
        double sum_q2=xy[0]*xy[0]+xy[1]*xy[1];
        if(sum_q2>radius*radius){
            ierr=1;
        }else{
            xyz[0]=2.0*xy[0]; xyz[1]=2.0*xy[1]; xyz[2]=1.0-sum_q2;
            vector_constant(xyz, 1.0/(radius*radius+sum_q2), xyz);
        }
    }
}

void vector_rotate(double r_v[3], int R[3][3], double v[3])
{
    double c_v[3];
    vector_copy(c_v, v);
    for(int i=0;i<3;i++){
        r_v[i]=0.0;
        for(int j=0;j<3;j++){
            r_v[i]+=(double)R[i][j]*c_v[j];
        }
    }
}

void image_pixels(const char* png_path, unsigned char *pixels, int numpx, int numpy)
{	
    png_structp png_ptr;  
    png_infop info_ptr;  
    FILE *fp=fopen(png_path, "wb");  
    if(!fp){
        printf("[ERROR] Unable to create %s using fopen.\n", png_path);
        exit(EXIT_FAILURE);
    }
    png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);  
    if(png_ptr==NULL){  
        printf("[ERROR] Unable to create %s using png_create_write_struct.\n", png_path);
        fclose(fp);
        exit(EXIT_FAILURE);
    }  
    info_ptr=png_create_info_struct(png_ptr);
    if(info_ptr==NULL){  
        printf("[ERROR] Unable to create %s using png_create_info_struct.\n", png_path);
        png_destroy_write_struct(&png_ptr, NULL);  
        exit(EXIT_FAILURE);
    }  
    png_init_io(png_ptr, fp);  
    png_set_IHDR(png_ptr, info_ptr, numpx, numpy, 8, 
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE); 
    png_colorp palette=(png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH*sizeof(png_color));
    if(!palette){
        printf("[ERROR] Unable to create palette using png_malloc.\n");
        fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        exit(EXIT_FAILURE);
    }
    png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);  
    png_write_info(png_ptr, info_ptr);  
    png_set_packing(png_ptr);
    png_bytepp rows=(png_bytepp)png_malloc(png_ptr, numpy*sizeof(png_bytep));
    for(int i=0;i<numpy;++i){
        rows[numpy-i-1]=(png_bytep)(pixels+(i)*numpx*3);
    }
    png_write_image(png_ptr, rows);  
    delete[] rows;
    png_write_end(png_ptr, info_ptr);  
    png_free(png_ptr, palette);  
    palette=NULL;  
    png_destroy_write_struct(&png_ptr, &info_ptr);  
    fclose(fp);   
}