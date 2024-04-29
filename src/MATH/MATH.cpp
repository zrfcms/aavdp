#include "MATH.h"

void image(const char* png_path, unsigned char *pixels, int width, int height)
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
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, 
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

void matrix_copy(double c_mat[3][3], double mat[3][3])
{
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            c_mat[i][j]=mat[i][j];
        }
    }
}

void matrix_constant(double k_mat[3][3], double k, double mat[3][3])
{
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            k_mat[i][j]=k*mat[j][i];
        }
    }
}

void matrix_transpose(double t_mat[3][3], double mat[3][3])
{
    double c_mat[3][3];
    matrix_copy(c_mat, mat);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            t_mat[i][j]=c_mat[j][i];
        }
    }
}

void matrix_multiply(double mat[3][3], double mat1[3][3], double mat2[3][3])
{
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            mat[i][j]=0.0;
            for(int k=0;k<3;k++){
                mat[i][j]+=mat1[i][k]*mat2[k][j];
            }
        }
    }
}

double vector_dot(double v1[3], double v2[3])
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double vector_length(double v[3])
{
    return sqrt(vector_dot(v, v));
}

void vector_normalize(double n_v[3], double v[3])
{
    double len_v=vector_length(v);
    if(0.0!=len_v){
        n_v[0]=v[0]/len_v; n_v[1]=v[1]/len_v; n_v[2]=v[2]/len_v;
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}

void vector_cross(double v[3], double v1[3], double v2[3])
{
    v[0]=v1[1]*v2[2]-v1[2]*v2[1];
    v[1]=v1[2]*v2[0]-v1[0]*v2[2];
    v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

void vector_constant(double k_v[3], double k, double v[3])
{
    k_v[0]=k*v[0]; k_v[1]=k*v[1]; k_v[2]=k*v[2];
}

void vector_plus(double v[3], double v1[3], double v2[3])
{
    v[0]=v1[0]+v2[0]; v[1]=v1[1]+v2[1]; v[2]=v1[2]+v2[2];
}

void vector_difference(double v[3], double v1[3], double v2[3])
{
    v[0]=v1[0]-v2[0]; v[1]=v1[1]-v2[1]; v[2]=v1[2]-v2[2];
}

void vector_multiply(double v[3], double v1[3], double v2[3])
{
    v[0]=v1[0]*v2[0]; v[1]=v1[1]*v2[1]; v[2]=v1[2]*v2[2];
}

void vector_rotate(double r_v[3], double R[3][3], double v[3])
{
    double c_v[3];
    vector_copy(c_v, v);
    for(int i=0;i<3;i++){
        r_v[i]=0.0;
        for(int j=0;j<3;j++){
            r_v[i]+=R[i][j]*c_v[j];
        }
    }
}

QUATERNION quate_conjg(QUATERNION q)
{
    QUATERNION qc={q.c1, -q.c2, -q.c3, -q.c4};
    return qc;
}

QUATERNION quate_multi(QUATERNION q1, QUATERNION q2)
{
    QUATERNION qm;
    double eps=1.0;
    qm.c1=q1.c1*q2.c1-q1.c2*q2.c2-q1.c3*q2.c3-q1.c4*q2.c4;
    qm.c2=q1.c1*q2.c2+q1.c2*q2.c1+eps*(q1.c3*q2.c4-q1.c4*q2.c3);
    qm.c3=q1.c1*q2.c3+q1.c3*q2.c1+eps*(q1.c4*q2.c2-q1.c2*q2.c4);
    qm.c4=q1.c1*q2.c4+q1.c4*q2.c1+eps*(q1.c2*q2.c3-q1.c3*q2.c2); 
    return qm;
}

void quate_rotate(double r_v[3], double v[3], QUATERNION q)
{
    QUATERNION qv={0.0, v[0], v[1], v[2]};
    QUATERNION r_qv=quate_multi(q, quate_multi(qv, quate_conjg(q)));
    r_v[0]=r_qv.c2; r_v[1]=r_qv.c3; r_v[2]=r_qv.c4; 
}

QUATERNION quate_convert(double v[3], double angle)
{
    QUATERNION q;
    if(fabs(angle)<1e-12){
        q.c1=1.0; q.c2=q.c3=q.c4=0.0;
    }else{
        double c=cos(0.5*angle);
        double s=sin(0.5*angle);
        q.c1=c; q.c2=s*v[0]; q.c3=s*v[1]; q.c4=s*v[2];
    }
    return q;
}