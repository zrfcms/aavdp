#include "./EBSD/EBSD_MASTER.h"
#include "./EBSD/EBSD_KKD.h"
#include "./EBSD/EBSD_CELL.h"
#include "./EBSD/EBSD_MC.h"
#include "./EBSD/EBSD_HDF5.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    // char nml_path[]="./test.nml";
    // char h5_path[]="./Ni-master.h5";
    // EBSD_MASTER master(nml_path);
    // master.compute_master_pattern(h5_path);

    EBSD_KKD kkd("./test.nml");
    //kkd.test_images("Ni-kkd-ref.h5");
    kkd.compute_master_pattern("Ni-kkd.h5");
}
