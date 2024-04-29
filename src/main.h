#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include "./MODEL/MODEL.h"
#include "./SED/SED.h"
#include "./KKD/KKD.h"
#include "./RDF/RDF.h"
#include "./SSF/SSF.h"
#include "./XRD/XRD.h"
#define PATH_CHAR_NUMBER 100
#define TYPE_INPUT_NUMBER 10

extern double be_double(char str[]);
extern bool is_parameter(char str[]);
extern void no_extension(char n_str[], char str[]);
extern void to_string(int num, char str[]);

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

void no_extension(char n_str[], char str[])
{
    char *ch=strrchr(str,'.');
    strcpy(n_str, str);
    if(ch!=nullptr){
        int pos=ch-str;
        n_str[pos]='\0';
    }
}

void to_string(int num, char str[])
{
    int i=0;
    if(num<0)
    {
        num=-num;
        str[i++]='-';
    } 
    do
    {
        str[i++]=num%10+48;//0-9: 48-57
        num /= 10;
    }while(num);
    str[i]='\0';
    int j=0;
    if(str[0]=='-')
    {
        j=1;
        ++i;
    } 
    for(;j<i/2;j++)
    {
        str[j]=str[j]+str[i-1-j];
        str[i-1-j]=str[j]-str[i-1-j];
        str[j]=str[j]-str[i-1-j];
    } 
}

#endif