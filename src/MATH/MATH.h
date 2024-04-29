#ifndef __AAVDP_MATH_H__
#define __AAVDP_MATH_H__
#include <cmath>
#include "../include/libpng/png.h"

#define KEY_CHAR_NUMBER 100
#define KEY_LINE_NUMBER 100

#define PI 3.141592653589793
#define DEG_TO_RAD 0.017453292519943295 //PI/180

extern void image(const char* png_path, unsigned char *pixels, int width, int height);

extern void matrix_copy(double c_mat[3][3], double mat[3][3]);
extern void matrix_constant(double k_mat[3][3], double k, double mat[3][3]);
extern void matrix_transpose(double t_mat[3][3], double mat[3][3]);
extern void matrix_multiply(double mat[3][3], double mat1[3][3], double mat2[3][3]);

template <typename T>
extern void vector_copy(T c_v[3], T v[3]);
extern double vector_dot(double v1[3], double v2[3]);
extern double vector_length(double v[3]);
extern void vector_normalize(double n_v[3], double v[3]);
extern void vector_cross(double v[3], double v1[3], double v2[3]);
extern void vector_constant(double k_v[3], double k, double v[3]);
extern void vector_plus(double v[3], double v1[3], double v2[3]);
extern void vector_difference(double v[3], double v1[3], double v2[3]);
extern void vector_multiply(double v[3], double v1[3], double v2[3]);
extern void vector_rotate(double r_v[3], double R[3][3], double v[3]);

template <typename T>
void vector_copy(T c_v[3], T v[3])
{
    c_v[0]=v[0]; c_v[1]=v[1]; c_v[2]=v[2];
}

struct QUATERNION
{
    double c1, c2, c3, c4;
};

extern QUATERNION quate_conjg(QUATERNION q);
extern QUATERNION quate_multi(QUATERNION q1, QUATERNION q2);
extern void quate_rotate(double r_v[3], double v[3], QUATERNION q);
extern QUATERNION quate_convert(double v[3], double angle);

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