#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include <unistd.h>
#include "./MODEL/MODEL.h"
#include "./XRD/XRD.h"
#include "./KED/KED.h"
#include "./DED/DED.h"
#include "./RDF/RDF.h"
#include "./SSF/SSF.h"

#define TYPE_INPUT_NUMBER 10
#define PATH_CHAR_NUMBER 100
#define EXT_CHAR_NUMBER 20

extern bool is_path_accessible(char path[]);
extern double be_double(char str[]);
extern bool is_parameter(char str[]);
extern void int_to_str(char str[], int num);
extern void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num);
extern void split_path(char name[], char ext[], char file_path[]);

bool is_path_accessible(char path[])
{
    if(access(path, F_OK)==-1){
        printf("[ERROR] Path %s does not exist", path);
        exit(1);
    }else if(access(path, R_OK)==-1){
        printf("[ERROR] Path %s is not readable", path);
        exit(1);
    }
    return true;
}

double be_double(char str[])
{
    char *err;
    double value=strtod(str, &err);
    if(strlen(err)!='\0')
    {
        printf("[ERROR] Input parameter %s should be float or int type", str);
        exit(1);
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

bool is_mode(char str[])
{
    if(0==strncmp(str, "--", 2)){
        return true;
    }
    return false;
}

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

void split_path(char name[], char ext[], char file_path[])
{
    char *ch=strrchr(file_path,'.');
    strcpy(name, file_path);
    strcpy(ext, ch);
    if(ch!=nullptr){
        int pos=ch-file_path;
        name[pos]='\0';
    }
}

#endif