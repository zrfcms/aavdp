#ifndef __AAVDP_MATH_H__
#define __AAVDP_MATH_H__
#include <cstdio>
#include <cmath>
#include <cstring>

#define PATH_CHAR_NUMBER 100
#define EXT_CHAR_NUMBER 20

#define ZERO_LIMIT 1.0e-4

#define PI 3.141592653589793
#define TWO_PI 6.283185307179586
#define FOUR_PI 12.566370614359172
#define DEG_TO_RAD 0.017453292519943295 //PI/180
#define RAD_TO_DEG 57.29577951308232 //180/PI

#define PI_SQRT_HALF 0.886226925452758 //sqrt(PI)/2
#define PI_HALF_SQRT 1.253314137315500 //sqrt(PI/2)
#define PI_INVERSE 0.318309886183791 //1/PI
#define SQRT_3_INVERSE 0.577350269190 //1/sqrt(3)

#define PI_PREE 0.7598356856515930 //3^(-1/4)
#define PI_PREF 1.3819765978853420 //sqrt(6/PI)
#define PI_PREA 0.5250375679043320 //3^(1/4)/sqrt(2*PI)
#define PI_PREG 1.55512030155621410 //2*sqrt(PI)/3^(3/4)
#define PI_PREC 0.9068996821171090 //PI/2*sqrt(3)
#define PI_PREB 1.0500751358086640 //3^(1/4)*sqrt(2/PI)
#define PI_PRED 2.0943951023931950 //2*PI/3
#define SQRT_3 1.732050807568877 //sqrt(3)

#define SQRT_HALF_3 0.866025403780 //sqrt(3)/2

#define ZERO 1e-12
#define ONE 0.99999
#define CONST_IN_AREA_INV 1.5273987E19 //1.0/(5.21E-21*(4*PI))

template <typename T>
extern T m_min(T a, T b);

template <typename T>
T m_min(T a, T b)
{
	return (a<b?a:b);
}

extern void int_to_str(char str[], int num);
extern void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num);
extern void split_path(char name[], char ext[], char file_path[]);

template <typename T>
extern void vector_copy(T c_v[3], T v[3]);
template <typename T>
extern void vector_constant(T k_v[3], T k, T v[3]);
template <typename T>
extern T vector_length(T v[3]);
template <typename T>
extern void vector_normalize(T n_v[3], T v[3]);
template <typename T>
extern T vector_dot(T v1[3], T v2[3]);
template <typename T>
extern void vector_cross(T v[3], T v1[3], T v2[3]);
template <typename T>
extern void vector_plus(T v[3], T v1[3], T v2[3]);
template <typename T>
extern void vector_difference(T v[3], T v1[3], T v2[3]);
template <typename T>
extern void vector_multiply(T v[3], T v1[3], T v2[3]);
template <typename T>
extern void vector_rotate(T r_v[3], T R[3][3], T v[3]);
extern void vector_rotate(double r_v[3], int R[3][3], double v[3]);
template <typename T>
extern void vector_transform(T t_v[3], T v[3], T R[3][3]);
template <typename T>
void vector_copy(T c_v[3], T v[3])
{
    c_v[0]=v[0]; c_v[1]=v[1]; c_v[2]=v[2];
}
template <typename T>
void vector_constant(T k_v[3], T k, T v[3])
{
    k_v[0]=k*v[0]; k_v[1]=k*v[1]; k_v[2]=k*v[2];
}
template <typename T>
T vector_length(T v[3])
{
    return sqrt(vector_dot(v, v));
}
template <typename T>
void vector_normalize(T n_v[3], T v[3])
{
    T len_v=vector_length(v);
    if(0.0!=len_v){
        n_v[0]=v[0]/len_v; n_v[1]=v[1]/len_v; n_v[2]=v[2]/len_v;
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}
template <typename T>
T vector_dot(T v1[3], T v2[3])
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}
template <typename T>
void vector_cross(T v[3], T v1[3], T v2[3])
{
    v[0]=v1[1]*v2[2]-v1[2]*v2[1];
    v[1]=v1[2]*v2[0]-v1[0]*v2[2];
    v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}
template <typename T>
void vector_plus(T v[3], T v1[3], T v2[3])
{
    v[0]=v1[0]+v2[0]; v[1]=v1[1]+v2[1]; v[2]=v1[2]+v2[2];
}
template <typename T>
void vector_difference(T v[3], T v1[3], T v2[3])
{
    v[0]=v1[0]-v2[0]; v[1]=v1[1]-v2[1]; v[2]=v1[2]-v2[2];
}
template <typename T>
void vector_multiply(T v[3], T v1[3], T v2[3])
{
    v[0]=v1[0]*v2[0]; v[1]=v1[1]*v2[1]; v[2]=v1[2]*v2[2];
}
template <typename T>
void vector_rotate(T r_v[3], T R[3][3], T v[3])
{
    T c_v[3];
    vector_copy(c_v, v);
    for(int i=0;i<3;i++){
        r_v[i]=0.0;
        for(int j=0;j<3;j++){
            r_v[i]+=R[i][j]*c_v[j];
        }
    }
}
template <typename T>
void vector_transform(T t_v[3], T v[3], T R[3][3])
{
    T c_v[3];
    vector_copy(c_v, v);
    for(int i=0;i<3;i++){
        t_v[i]=0.0;
        for(int j=0;j<3;j++){
            t_v[i]+=R[j][i]*c_v[j];
        }
    }
}

template <typename T>
extern void matrix_copy(T c_mat[3][3], T mat[3][3]);
template <typename T>
extern void matrix_constant(T k_mat[3][3], T k, T mat[3][3]);
template <typename T>
extern void matrix_transpose(T t_mat[3][3], T mat[3][3]);
template <typename T>
extern void matrix_multiply(T mat[3][3], T mat1[3][3], T mat2[3][3]);
template <typename T>
void matrix_copy(T c_mat[3][3], T mat[3][3])
{
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            c_mat[i][j]=mat[i][j];
        }
    }
}
template <typename T>
void matrix_constant(T k_mat[3][3], T k, T mat[3][3])
{
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            k_mat[i][j]=k*mat[i][j];
        }
    }
}
template <typename T>
void matrix_transpose(T t_mat[3][3], T mat[3][3])
{
    T c_mat[3][3];
    matrix_copy(c_mat, mat);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            t_mat[i][j]=c_mat[j][i];
        }
    }
}
template <typename T>
void matrix_multiply(T mat[3][3], T mat1[3][3], T mat2[3][3])
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

extern void compute_square_Lambert(double xy[2], int &ierr, double xyz[3]);
extern void compute_sphere_from_square_Lambert(double xyz[3], int &ierr, double xy[2]);
extern  int get_sextant(double x, double y);
extern void compute_hexagonal_Lambert(double xy[2], int &ierr, double xyz[3]);
extern void compute_sphere_from_hexagonal_Lambert(double xyz[3], int &ierr, double xy[2]);
extern void compute_sphere_from_stereographic_projection(double xyz[3], int &ierr, double xy[2], double radius=1.0);

template <typename T>
extern void mallocate(T **data, int size);
template <typename T>
extern void callocate(T **data, int size, T value);
template <typename T>
extern void deallocate(T *data);
template <typename T>
extern void mallocate_2d(T ***data, int nrow, int ncol);
template <typename T>
extern void callocate_2d(T ***data, int nrow, int ncol, T value);
template <typename T>
extern void deallocate_2d(T **data, int nrow);
template <typename T>
extern void mallocate_3d(T ****data, int nrow, int ncol, int nlayer);
template <typename T>
extern void callocate_3d(T ****data, int nrow, int ncol, int nlayer, T value);
template <typename T>
extern void deallocate_3d(T ***data, int nrow, int ncol);
template <typename T>
extern void mallocate_4d(T *****data, int nrow, int ncol, int nlayer, int nblock);
template <typename T>
extern void callocate_4d(T *****data, int nrow, int ncol, int nlayer, int nblock, T value);
template <typename T>
extern void deallocate_4d(T ****data, int nrow, int ncol, int nlayer);
template <typename T>
extern void reshape_2d(T ***data, T *wdata, int nrow, int ncol);
template <typename T>
extern void reshape_3d(T ****data, T *wdata, int nrow, int ncol, int nlayer);
template <typename T>
extern void reshape_4d(T *****data, T *wdata, int nrow, int ncol, int nlayer, int nblock);
template <typename T>
extern void unreshape_2d(T **wdata, T **data, int nrow, int ncol);
template <typename T>
extern void unreshape_3d(T **wdata, T ***data, int nrow, int ncol, int nlayer);
template <typename T>
extern void unreshape_4d(T **wdata, T ****data, int nrow, int ncol, int nlayer, int nblock);
template <typename T>
void mallocate(T **data, int size)
{
    (*data)=new T[size];
}
template <typename T>
void callocate(T **data, int size, T value)
{
    mallocate(data, size);
    for(int i=0;i<size;i++){
        (*data)[i]=value;
    }
}
template <typename T>
void deallocate(T *data)
{
    if(data==nullptr) return;
    delete [] data;
}
template <typename T>
void mallocate_2d(T ***data, int nrow, int ncol)
{
    (*data)=new T*[nrow];
    for(int i=0;i<nrow;i++){
        (*data)[i]=new T[ncol];
    }
}
template <typename T>
void callocate_2d(T ***data, int nrow, int ncol, T value)
{
    mallocate_2d(data, nrow, ncol);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            (*data)[i][j]=value;
        }
    }
}
template <typename T>
void deallocate_2d(T **data, int nrow)
{
    if(data==nullptr) return;
    for(int i=0;i<nrow;i++){
        delete [] data[i];
    }
    delete [] data;
}
template <typename T>
void mallocate_3d(T ****data, int nrow, int ncol, int nlayer)
{
    (*data)=new T**[nrow];
    for(int i=0;i<nrow;i++){
        (*data)[i]=new T*[ncol];
    }
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            (*data)[i][j]=new T[nlayer];
        }
    }
}
template <typename T>
void callocate_3d(T ****data, int nrow, int ncol, int nlayer, T value)
{
    mallocate_3d(data, nrow, ncol, nlayer);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                (*data)[i][j][k]=value;
            }
        }
    }
}
template <typename T>
void deallocate_3d(T ***data, int nrow, int ncol)
{
    if(data==nullptr) return;
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            delete [] data[i][j];
        }
    }
    for(int i=0;i<nrow;i++){
        delete [] data[i];
    }
    delete [] data;
}
template <typename T>
void mallocate_4d(T *****data, int nrow, int ncol, int nlayer, int nblock)
{
    (*data)=new T***[nrow];
    for(int i=0;i<nrow;i++){
        (*data)[i]=new T**[ncol];
    }
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            (*data)[i][j]=new T*[nlayer];
        }
    }
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                (*data)[i][j][k]=new T[nblock];
            }
        }
    }
}
template <typename T>
void callocate_4d(T *****data, int nrow, int ncol, int nlayer, int nblock, T value)
{
    mallocate_4d(data, nrow, ncol, nlayer, nblock);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                for(int n=0;n<nblock;n++){
                    (*data)[i][j][k][n]=value;
                }
            }
        }
    }
}
template <typename T>
void deallocate_4d(T ****data, int nrow, int ncol, int nlayer)
{
    if(data==nullptr) return;
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                delete [] data[i][j][k];
            }
        }
    }
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            delete [] data[i][j];
        }
    }
    for(int i=0;i<nrow;i++){
        delete [] data[i];
    }
    delete [] data;
}
template <typename T>
void reshape_2d(T ***data, T *wdata, int nrow, int ncol)
{
    mallocate_2d(data, nrow, ncol);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            (*data)[i][j]=wdata[i*ncol+j];
        }
    }
}
template <typename T>
void unreshape_2d(T **wdata, T **data, int nrow, int ncol)
{
    mallocate(wdata, nrow*ncol);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            (*wdata)[i*ncol+j]=data[i][j];
        }
    }
}
template <typename T>
void reshape_3d(T ****data, T *wdata, int nrow, int ncol, int nlayer)
{
    mallocate_3d(data, nrow, ncol, nlayer);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                (*data)[i][j][k]=wdata[i*ncol*nlayer+j*nlayer+k];
            }
        }
    }
}
template <typename T>
void unreshape_3d(T **wdata, T ***data, int nrow, int ncol, int nlayer)
{
    mallocate(wdata, nrow*ncol*nlayer);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                (*wdata)[i*ncol*nlayer+j*nlayer+k]=data[i][j][k];
            }
        }
    }
}
template <typename T>
void reshape_4d(T *****data, T *wdata, int nrow, int ncol, int nlayer, int nblock)
{
    mallocate_4d(data, nrow, ncol, nlayer, nblock);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                for(int n=0;n<nblock;n++){
                    (*data)[i][j][k][n]=wdata[i*ncol*nlayer*nblock+j*nlayer*nblock+k*nblock+n];
                }
            }
        }
    }
}
template <typename T>
void unreshape_4d(T **wdata, T ****data, int nrow, int ncol, int nlayer, int nblock)
{
    mallocate(wdata, nrow*ncol*nlayer*nblock);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<nlayer;k++){
                for(int n=0;n<nblock;n++){
                    (*wdata)[i*ncol*nlayer*nblock+j*nlayer*nblock+k*nblock+n]=data[i][j][k][n];
                }
            }
        }
    }
}
#endif