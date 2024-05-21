#ifndef _AAVDP_HDF5_H_
#define _AAVDP_HDF5_H_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "../../include/hdf5/hdf5.h"
#include "../MATH/MATH.h"
using namespace std;

#define HDF5_STRING_LENGTH 100

class HDF5
{
public:
    void open(const char *file_path);
    void close();
    void write_group(const char *name);
    void write(const char *name, int data);
    void write(const char *name, double data);
    void write(const char *name, char *data);
    void write_array(const char *name, int *data, size_t size);
    void write_array(const char *name, double *data, size_t size);
    void write_array_2d(const char *name, int **data, size_t size1, size_t size2);
    void write_array_2d(const char *name, double **data, size_t size1, size_t size2);
    void write_array_3d(const char *name, int ***data, size_t size1, size_t size2, size_t size3);
    void write_array_3d(const char *name, double ***data, size_t size1, size_t size2, size_t size3);
    void write_hyper_array_3d(const char *name, double ***data, size_t *mem_dims, size_t *mem_offsets, size_t size1, size_t size2, size_t size3);
    void write_array_4d(const char *name, int ****data, size_t size1, size_t size2, size_t size3, size_t size4);
    void write_array_4d(const char *name, double ****data, size_t size1, size_t size2, size_t size3, size_t size4);
    void write_hyper_array_4d(const char *name, double ****data, size_t *mem_dims, size_t *mem_offsets, size_t size1, size_t size2, size_t size3, size_t size4);
    void read(const char *name, int &data);
    void read(const char *name, double &data);
    void read(const char *name, char data[]);
    void read_array(const char *name, int **data, size_t &size);
    void read_array(const char *name, double **data, size_t &size);
    void read_array_2d(const char *name, int ***data, size_t &size1, size_t &size2);
    void read_array_2d(const char *name, double ***data, size_t &size1, size_t &size2);
    void read_array_3d(const char *name, int ****data, size_t &size1, size_t &size2, size_t &size3);
    void read_array_3d(const char *name, double ****data, size_t &size1, size_t &size2, size_t &size3);
    void read_array_4d(const char *name, int *****data, size_t &size1, size_t &size2, size_t &size3, size_t &size4);
    void read_array_4d(const char *name, double *****data, size_t &size1, size_t &size2, size_t &size3, size_t &size4);
private:
    hid_t  f_id;
    herr_t status;
    hid_t  dataspace_id, dataset_id;
    hid_t  type_id;
    void   create_dataset(const char *name, size_t ndim=0, size_t size1=1, size_t size2=1, size_t size3=1, size_t size4=1);
    void   open_dataset(const char *name, hsize_t *dims, unsigned ndim=0);
    void   close_dataset();
};

#endif