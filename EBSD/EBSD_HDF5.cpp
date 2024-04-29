#include "EBSD_HDF5.h"

void EBSD_HDF5::open(const char* file_path)
{
    status=H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    bool existing_flag=false;
    if(FILE *fp=fopen(file_path, "r")){
        fclose(fp);
        existing_flag=true;
    }
    if(!existing_flag){
        f_id=H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if(f_id<0){
            printf("[ERROR] Unable to create HDF5 file %s.\n", file_path);
            exit(EXIT_FAILURE);
        }
    }else{
        f_id=H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);
        if(f_id<0){
            printf("[ERROR] Unable to open HDF5 file %s.\n", file_path);
            exit(EXIT_FAILURE);
        }
    }
}

void EBSD_HDF5::close()
{
    status=H5Fclose(f_id);
    status=H5Eset_auto(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr);
}

void EBSD_HDF5::write_group(const char *name)
{
    hid_t g_id;
    status=H5Lget_info(f_id, name, NULL, H5P_DEFAULT);
    if(0!=status){
        g_id=H5Gcreate(f_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(g_id<0){
            printf("[ERROR] Unable to create HDF5 group %s.\n", name);
            exit(EXIT_FAILURE);
        }
        status=H5Gclose(g_id);
    }
}

void EBSD_HDF5::write(const char *name, int data)
{
    type_id=H5T_NATIVE_INT;
    create_dataset(name);
    status=H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    close_dataset();
}

void EBSD_HDF5::write(const char *name, double data)
{
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    close_dataset();
}

void EBSD_HDF5::write(const char *name, char *data)
{
    type_id=H5Tcopy(H5T_C_S1);
    status=H5Tset_size(type_id, H5T_VARIABLE);
    create_dataset(name, 0, strlen(data));
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    close_dataset();
}

void EBSD_HDF5::write_array(const char *name, int *data, size_t size)
{
    type_id=H5T_NATIVE_INT;
    create_dataset(name, 1, size);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    close_dataset();
}

void EBSD_HDF5::write_array(const char *name, double *data, size_t size)
{
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name, 1, size);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    close_dataset();
}

void EBSD_HDF5::write_array_2d(const char *name, int **data, size_t size1, size_t size2)
{
    int *wdata;
    reshape2d(&wdata, data, size1, size2);
    type_id=H5T_NATIVE_INT;
    create_dataset(name, 2, size1, size2);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_array_2d(const char *name, double **data, size_t size1, size_t size2)
{
    double *wdata;
    reshape2d(&wdata, data, size1, size2);
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name, 2, size1, size2);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_array_3d(const char *name, int ***data, size_t size1, size_t size2, size_t size3)
{
    int *wdata;
    reshape3d(&wdata, data, size1, size2, size3);
    type_id=H5T_NATIVE_INT;
    create_dataset(name, 3, size1, size2, size3);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_array_3d(const char *name, double ***data, size_t size1, size_t size2, size_t size3)
{
    double *wdata;
    reshape3d(&wdata, data, size1, size2, size3);
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name, 3, size1, size2, size3);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_hyper_array_3d(const char *name, double ***data, size_t *mem_dims, size_t *mem_offsets, size_t size1, size_t size2, size_t size3)
{
    double *wdata;
    reshape3d(&wdata, data, mem_dims[0], mem_dims[1], mem_dims[2]);
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name, 3, size1, size2, size3);
    hsize_t mdims[4]={mem_dims[0], mem_dims[1], mem_dims[2]};
    hsize_t moffsets[4]={mem_offsets[0], mem_offsets[1], mem_offsets[2]};
    hid_t memspace_id=H5Screate_simple(3, mdims, NULL);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, moffsets, NULL, mdims, NULL);
    status=H5Dwrite(dataset_id, type_id, memspace_id, dataspace_id, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_array_4d(const char *name, int ****data, size_t size1, size_t size2, size_t size3, size_t size4)
{
    int *wdata;
    reshape4d(&wdata, data, size1, size2, size3, size4);
    type_id=H5T_NATIVE_INT;
    create_dataset(name, 4, size1, size2, size3, size4);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_array_4d(const char *name, double ****data, size_t size1, size_t size2, size_t size3, size_t size4)
{
    double *wdata;
    reshape4d(&wdata, data, size1, size2, size3, size4);
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name, 4, size1, size2, size3, size4);
    status=H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::write_hyper_array_4d(const char *name, double ****data, size_t *mem_dims, size_t *mem_offsets, size_t size1, size_t size2, size_t size3, size_t size4)
{
    double *wdata;
    reshape4d(&wdata, data, mem_dims[0], mem_dims[1], mem_dims[2], mem_dims[3]);
    type_id=H5T_NATIVE_DOUBLE;
    create_dataset(name, 4, size1, size2, size3, size4);
    hsize_t mdims[4]={mem_dims[0], mem_dims[1], mem_dims[2], mem_dims[3]};
    hsize_t moffsets[4]={mem_offsets[0], mem_offsets[1], mem_offsets[2], mem_offsets[3]};
    hid_t memspace_id=H5Screate_simple(4, mdims, NULL);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, moffsets, NULL, mdims, NULL);
    status=H5Dwrite(dataset_id, type_id, memspace_id, dataspace_id, H5P_DEFAULT, wdata);
    deallocate(wdata);
    close_dataset();
}

void EBSD_HDF5::read(const char *name, int &data)
{
    hsize_t dims[1]={0};
    open_dataset(name, dims);
    type_id=H5T_NATIVE_INT;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    close_dataset();
}

void EBSD_HDF5::read(const char *name, double &data)
{
    hsize_t dims[1]={0};
    open_dataset(name, dims);
    type_id=H5T_NATIVE_DOUBLE;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    close_dataset();
}

void EBSD_HDF5::read(const char *name, char *data)
{
    hsize_t dims[1]={0};
    open_dataset(name, dims);
    type_id=H5Dget_type(dataset_id);
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    close_dataset();
}

void EBSD_HDF5::read_array(const char *name, int **data, size_t &size)
{
    unsigned ndim=1;
    hsize_t dims[1]={0}; 
    open_dataset(name, dims, ndim); 
    size=dims[0];
    mallocate(data, size);
    type_id=H5T_NATIVE_INT;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (*data));
    close_dataset();
}

void EBSD_HDF5::read_array(const char *name, double **data, size_t &size)
{
    unsigned ndim=1;
    hsize_t dims[1]={0}; 
    open_dataset(name, dims, ndim); 
    size=dims[0];
    mallocate(data, size);
    type_id=H5T_NATIVE_DOUBLE;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (*data));
    close_dataset();
}

void EBSD_HDF5::read_array_2d(const char *name, int ***data, size_t &size1, size_t &size2)
{
    unsigned ndim=2;
    hsize_t dims[2]={0}; 
    open_dataset(name, dims, ndim); 
    size1=dims[0]; size2=dims[1];
    int *wdata=nullptr; 
    mallocate(&wdata, size1*size2);
    type_id=H5T_NATIVE_INT;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    close_dataset();
    reshape_2d(data, wdata, size1, size2);
}

void EBSD_HDF5::read_array_2d(const char *name, double ***data, size_t &size1, size_t &size2)
{
    unsigned ndim=2;
    hsize_t dims[2]={0}; 
    open_dataset(name, dims, ndim); 
    size1=dims[0]; size2=dims[1];
    double *wdata=nullptr; 
    mallocate(&wdata, size1*size2);
    type_id=H5T_NATIVE_DOUBLE;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    close_dataset();
    reshape_2d(data, wdata, size1, size2);
}

void EBSD_HDF5::read_array_3d(const char *name, int ****data, size_t &size1, size_t &size2, size_t &size3)
{
    unsigned ndim=3;
    hsize_t dims[3]={0}; 
    open_dataset(name, dims, ndim); 
    size1=dims[0]; size2=dims[1]; size3=dims[2];
    int *wdata=nullptr; 
    mallocate(&wdata, size1*size2*size3);
    type_id=H5T_NATIVE_INT;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    close_dataset();
    reshape_3d(data, wdata, size1, size2, size3);
}

void EBSD_HDF5::read_array_3d(const char *name, double ****data, size_t &size1, size_t &size2, size_t &size3)
{
    unsigned ndim=3;
    hsize_t dims[3]={0}; 
    open_dataset(name, dims, ndim); 
    size1=dims[0]; size2=dims[1]; size3=dims[2];
    double *wdata=nullptr; 
    mallocate(&wdata, size1*size2*size3);
    type_id=H5T_NATIVE_DOUBLE;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    close_dataset();
    reshape_3d(data, wdata, size1, size2, size3);
}

void EBSD_HDF5::read_array_4d(const char *name, int *****data, size_t &size1, size_t &size2, size_t &size3, size_t &size4)
{
    unsigned ndim=4;
    hsize_t dims[4]={0}; 
    open_dataset(name, dims, ndim); 
    size1=dims[0]; size2=dims[1]; size3=dims[2]; size4=dims[3];
    int *wdata=nullptr; 
    mallocate(&wdata, size1*size2*size3*size4);
    type_id=H5T_NATIVE_INT;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    close_dataset();
    reshape_4d(data, wdata, size1, size2, size3, size4);
}

void EBSD_HDF5::read_array_4d(const char *name, double *****data, size_t &size1, size_t &size2, size_t &size3, size_t &size4)
{
    unsigned ndim=4;
    hsize_t dims[4]={0}; 
    open_dataset(name, dims, ndim); 
    size1=dims[0]; size2=dims[1]; size3=dims[2]; size4=dims[3];
    double *wdata=nullptr; 
    mallocate(&wdata, size1*size2*size3*size4);
    type_id=H5T_NATIVE_DOUBLE;
    status=H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    close_dataset();
    reshape_4d(data, wdata, size1, size2, size3, size4);
}

void EBSD_HDF5::create_dataset(const char *name, size_t ndim, size_t size1, size_t size2, size_t size3, size_t size4)
{
    status=H5Lget_info(f_id, name, NULL, H5P_DEFAULT); 
    if(0==status){
        hsize_t sizes[4]={size1, size2, size3, size4};
        hsize_t dims[ndim];
        bool flag=false;
        open_dataset(name, dims, ndim);
        for(int i=0;i<ndim;i++){
            if(dims[i]!=sizes[i]){
                flag=true;
            }
        }
        if(flag){
            printf("[ERROR] Unable to create HDF5 dataset %s due to the existing dataset with different dimensions.\n", name);
            exit(EXIT_FAILURE);
        }
    }else{
        size_t rank;
        hsize_t *dims;
        bool is_string_flag=false;
        switch(ndim)
        {
        case 0:
            rank=1; 
            mallocate(&dims, 1); dims[0]=size1;
            if(size1>1) is_string_flag=true;
            break;
        case 1:
            rank=1;
            mallocate(&dims, 1); dims[0]=size1;
            break;
        case 2:
            rank=2;
            mallocate(&dims, 2); dims[0]=size1; dims[1]=size2;
            break;
        case 3:
            rank=3;
            mallocate(&dims, 3); dims[0]=size1; dims[1]=size2; dims[2]=size3;
            break;
        case 4:
            rank=4;
            mallocate(&dims, 4); dims[0]=size1; dims[1]=size2; dims[2]=size3; dims[3]=size4;
            break;
        default:
            printf("[ERROR] Unable to store dataset of %d dimensions as HDF5 format.", ndim);
            exit(EXIT_FAILURE);
        }
        if(is_string_flag){
            dataspace_id=H5Screate(H5S_SCALAR);
        }else{
            dataspace_id=H5Screate_simple(rank, dims, NULL);
        }
        if(dataspace_id<0){
            printf("[ERROR] Unable to create HDF5 dataspace %s.\n", name);
            exit(EXIT_FAILURE);
        }
        dataset_id=H5Dcreate(f_id, name, type_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(dataset_id<0){
            printf("[ERROR] Unable to create HDF5 dataset %s.\n", name);
            exit(EXIT_FAILURE);
        }
    }
}

void EBSD_HDF5::open_dataset(const char *name, hsize_t *dims, unsigned ndim)
{
    dataset_id=H5Dopen(f_id, name, H5P_DEFAULT);
    dataspace_id=H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    unsigned rank=H5Sget_simple_extent_ndims(dataspace_id);
    if((rank==ndim)||((1==rank)&&(0==ndim))){
        return;
    }else{
        printf("[ERROR] Inconsistency between ndim %d preseted and rank %d deduced from dataset %s.\n", ndim, rank, name);
        exit(EXIT_FAILURE);
    }
}

void EBSD_HDF5::close_dataset()
{
    status=H5Dclose(dataset_id);
    status=H5Sclose(dataspace_id);
}


