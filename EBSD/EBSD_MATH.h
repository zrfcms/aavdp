#ifndef __EBSD_MATH_H__
#define __EBSD_MATH_H__
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "EBSD_CONST.h"
#include "../include/libpng/png.h"
using namespace std;

template <typename T>
extern void mallocate(T **data, int size);
template <typename T>
extern void callocate(T **data, int size, T value);
template <typename T>
extern void deallocate(T *data);
template <typename T>
extern void mallocate_2d(T ***data, int size1, int size2);
template <typename T>
extern void callocate_2d(T ***data, int size1, int size2, T value);
template <typename T>
extern void deallocate_2d(T **data, int size1);
template <typename T>
extern void mallocate_3d(T ****data, int size1, int size2, int size3);
template <typename T>
extern void callocate_3d(T ****data, int size1, int size2, int size3, T value);
template <typename T>
extern void deallocate_3d(T ***data, int size1, int size2);
template <typename T>
extern void mallocate_4d(T *****data, int size1, int size2, int size3, int size4);
template <typename T>
extern void callocate_4d(T *****data, int size1, int size2, int size3, int size4, T value);
template <typename T>
extern void deallocate_4d(T ****data, int size1, int size2, int size3);
template <typename T>
extern void reshape_2d(T ***data, T *wdata, size_t size1, size_t size2);
template <typename T>
extern void reshape_3d(T ****data, T *wdata, size_t size1, size_t size2, size_t size3);
template <typename T>
extern void reshape_4d(T *****data, T *wdata, size_t size1, size_t size2, size_t size3, size_t size4);
template <typename T>
extern void reshape2d(T **wdata, T **data, size_t size1, size_t size2);
template <typename T>
extern void reshape3d(T **wdata, T ***data, size_t size1, size_t size2, size_t size3);
template <typename T>
extern void reshape4d(T **wdata, T ****data, size_t size1, size_t size2, size_t size3, size_t size4);

extern bool is_in_vector(const int value, const int vec[], const int num);
template <typename T>
extern void normalize_vector(T n_v[], T v[], int num=3);

extern void compute_square_Lambert(double xy[2], int &ierr, double xyz[3]);
extern void compute_hexagonal_Lambert(double xy[2], int &ierr, double xyz[3]);
extern void compute_Lambert_interpolation(double xyz[3], int nump, bool hexagonal_flag, int &ix, int &iy, int &ixp, int &iyp, double &dx, double &dy, double &dxm, double &dym);

extern int get_sextant(double x, double y);
extern void compute_sphere_from_square_Lambert(double xyz[3], int &ierr, double xy[2]);
extern void compute_sphere_from_hexagonal_Lambert(double xyz[3], int &ierr, double xy[2]);
extern void compute_sphere_from_stereographic_projection(double xyz[3], int &ierr, double xy[2], double radius=1.0);

extern void create_image(const char* png_path, unsigned char *pixels, int width, int height);

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
void mallocate_2d(T ***data, int size1, int size2)
{
    (*data)=new T*[size1];
    for(int i=0;i<size1;i++){
        (*data)[i]=new T[size2];
    }
}

template <typename T>
void callocate_2d(T ***data, int size1, int size2, T value)
{
    mallocate_2d(data, size1, size2);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            (*data)[i][j]=value;
        }
    }
}

template <typename T>
void deallocate_2d(T **data, int size1)
{
    if(data==nullptr) return;
    for(int i=0;i<size1;i++){
        delete [] data[i];
    }
    delete [] data;
}

template <typename T>
void mallocate_3d(T ****data, int size1, int size2, int size3)
{
    (*data)=new T**[size1];
    for(int i=0;i<size1;i++){
        (*data)[i]=new T*[size2];
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            (*data)[i][j]=new T[size3];
        }
    }
}

template <typename T>
void callocate_3d(T ****data, int size1, int size2, int size3, T value)
{
    mallocate_3d(data, size1, size2, size3);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                (*data)[i][j][k]=value;
            }
        }
    }
}

template <typename T>
void deallocate_3d(T ***data, int size1, int size2)
{
    if(data==nullptr) return;
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            delete [] data[i][j];
        }
    }
    for(int i=0;i<size1;i++){
        delete [] data[i];
    }
    delete [] data;
}

template <typename T>
void mallocate_4d(T *****data, int size1, int size2, int size3, int size4)
{
    (*data)=new T***[size1];
    for(int i=0;i<size1;i++){
        (*data)[i]=new T**[size2];
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            (*data)[i][j]=new T*[size3];
        }
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                (*data)[i][j][k]=new T[size4];
            }
        }
    }
}

template <typename T>
void callocate_4d(T *****data, int size1, int size2, int size3, int size4, T value)
{
    mallocate_4d(data, size1, size2, size3, size4);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                for(int n=0;n<size4;n++){
                    (*data)[i][j][k][n]=value;
                }
            }
        }
    }
}

template <typename T>
void deallocate_4d(T ****data, int size1, int size2, int size3)
{
    if(data==nullptr) return;
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                delete [] data[i][j][k];
            }
        }
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            delete [] data[i][j];
        }
    }
    for(int i=0;i<size1;i++){
        delete [] data[i];
    }
    delete [] data;
}

template <typename T>
void reshape_2d(T ***data, T *wdata, size_t size1, size_t size2)
{
    mallocate_2d(data, size1, size2);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            (*data)[i][j]=wdata[i*size2+j];
        }
    }
}

template <typename T>
void reshape2d(T **wdata, T **data, size_t size1, size_t size2)
{
    mallocate(wdata, size1*size2);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            (*wdata)[i*size2+j]=data[i][j];
        }
    }
}

template <typename T>
void reshape_3d(T ****data, T *wdata, size_t size1, size_t size2, size_t size3)
{
    mallocate_3d(data, size1, size2, size3);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                (*data)[i][j][k]=wdata[i*size2*size3+j*size3+k];
            }
        }
    }
}

template <typename T>
void reshape3d(T **wdata, T ***data, size_t size1, size_t size2, size_t size3)
{
    mallocate(wdata, size1*size2*size3);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                (*wdata)[i*size2*size3+j*size3+k]=data[i][j][k];
            }
        }
    }
}

template <typename T>
void reshape_4d(T *****data, T *wdata, size_t size1, size_t size2, size_t size3, size_t size4)
{
    mallocate_4d(data, size1, size2, size3, size4);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                for(int n=0;n<size4;n++){
                    (*data)[i][j][k][n]=wdata[i*size2*size3*size4+j*size3*size4+k*size4+n];
                }
            }
        }
    }
}

template <typename T>
void reshape4d(T **wdata, T ****data, size_t size1, size_t size2, size_t size3, size_t size4)
{
    mallocate(wdata, size1*size2*size3*size4);
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                for(int n=0;n<size4;n++){
                    (*wdata)[i*size2*size3*size4+j*size3*size4+k*size4+n]=data[i][j][k][n];
                }
            }
        }
    }
}


#endif