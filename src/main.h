#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include "./MODEL/MODEL.h"
#include "./XRD/XRD.h"
#include "./SED/SED.h"
#include "./KKD/KKD.h"
#include "./DKD/DKD_MC.h"
#include "./DKD/DKD.h"
#include "./RDF/RDF.h"
#include "./SSF/SSF.h"

extern double be_double(char str[]);
extern bool is_parameter(char str[]);

double be_double(char str[])
{
    char *err;
    double value=strtod(str, &err);
    if(strlen(err)!='\0')
    {
        printf("[ERROR] Input parameter %s should be float or int type", str);
        exit(EXIT_FAILURE);
    }
    return value;
}

bool is_parameter(char str[])
{
    if(0==strncmp(str, "-", 1)){
        return false;
    }
    return true;
}

#endif