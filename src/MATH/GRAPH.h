#ifndef __AAVDP_GRAPH_H__
#define __AAVDP_GRAPH_H__
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "MATH.h"
#include "../include/png.h"

extern void image_pixels(const char* png_path, unsigned char *pixels, int numpx, int numpy);
template <typename T>
extern void image_array(const char* png_path, T **arr, int nrow, int ncol, double width, double height, int resolution, bool is_black_background=true);
template <typename T>
void image_array(const char* png_path, T **arr, int nrow, int ncol, double width, double height, int resolution, bool is_black_background)
{
    int numpx=width*resolution, numpy=height*resolution;
    int multiplier=m_min(numpx/ncol, numpy/nrow);
    numpx=ncol*multiplier; numpy=nrow*multiplier;
    double **warr; callocate_2d(&warr, numpy, numpx, 0.0);
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            for(int k=0;k<multiplier;k++){
                for(int n=0;n<multiplier;n++){
                    warr[i*multiplier+k][j*multiplier+n]=arr[i][j];
                }
            }
        }
    }
    int nump=numpx*numpy;
    double *wdata; unreshape_2d(&wdata, warr, numpy, numpx);
    double wmax=wdata[0], wmin=wdata[0];
    for(int i=0;i<nump;i++){
        if(wmax<wdata[i]) wmax=wdata[i];
        if(wmin>wdata[i]) wmin=wdata[i];
    }
    double wdiff=wmax-wmin;
    double wref=is_black_background?wmin:wmax;

    unsigned char *pixels;
    mallocate(&pixels, 3*nump);
    for(int i=0;i<nump;i++){
        pixels[i*3]=pixels[i*3+1]=pixels[i*3+2]=round(fabs(wref-wdata[i])/wdiff*255.0);
    }
    image_pixels(png_path, pixels, numpx, numpy);
}

struct POINT
{
    double x=0.0, y=0.0;
    double black=1.0;
    char   style='c';
    double psize=5.0;
    bool   is_in_area_checked=true;
};

struct AXIS
{
    int    n_major_tick;
    double *major_ticks;
    double major_tick_psize=20.0;
    bool   is_tick_in=true;
    bool   is_visible=true;
    double black=1.0;
    double line_pwidth=5.0;
    bool   is_in_area_checked=false;
};

struct AREA
{
    double border_width_factor=0.1;
    double border_height_factor=0.1;
    int    pwidth, pheight;
    int    start_px, start_py;
    int    end_px, end_py;
    double spacing_x, spacing_y;
    double start_x, end_x;
    double start_y, end_y;
    double ratio_x, ratio_y;
};

struct IMAGE
{
    int    pwidth, pheight;
    double black=0.0;
    double **pixels=nullptr;
};

class GRAPH
{
public:
    char   style[10];
    IMAGE  image;
    AREA   area;
    AXIS   left, right, top, down;
    int    nump=0;
    POINT  *points;
    GRAPH(double width, double height, int resolution);
    ~GRAPH();
    void set_xlim(double xmin, double xmax);
    void set_ylim(double ymin, double ymax);
    void set_xticks(double *major_ticks, double n_major_tick);
    void set_yticks(double *major_ticks, double n_major_tick);
    void set_tick_in(bool is_tick_in);
    void set_top_visible(bool is_visible);
    void set_right_visible(bool is_visible);
    void scatter(double *x, double *y, double *value, int num, char marker_style='c', double marker_size=20.0);
    void line(double *x, double *y, int num, double line_width=5.0);
    void hist(double *x, double *y, int num, double line_width=5.0);
    void draw(const char *png_path);
private:
    void get_pixel(int &i, int &j, double x, double y);
    void draw_marker(POINT *point);
    void draw_line(POINT *start_point, POINT *end_point);
    void draw_xaxis(AXIS *axis, double x1, double x2, double y);
    void draw_yaxis(AXIS *axis, double y1, double y2, double x);
    // inline void rgb_to_hsl(int hsl[3], int rgb[3]){
    //     double r=rgb[0]/255.0;
    //     double g=rgb[1]/255.0;
    //     double b=rgb[2]/255.0;
    //     double max_rgb=max3(r, g, b);
    //     double min_rgb=min3(r, g, b);

    //     double h, s, l;
    //     int    hsl[3];
    //     l=(max_rgb+min_rgb)/2.0;
    //     if(l<0.5){
    //         s=(max_rgb-min_rgb)/(max_rgb+min_rgb);
    //     }else{
    //         s=(max_rgb-min_rgb)/(2.0-max_rgb-min_rgb);
    //     }
    //     if(r==max_rgb){
    //         h=60.0*(g-b)/(max_rgb-min_rgb);
    //     }else if(g==max_rgb){
    //         h=60.0*((b-r)/(max_rgb-min_rgb)+2.0);
    //     }else{
    //         h=60.0*((r-g)/(max_rgb-min_rgb)+4.0);
    //     }
    //     hsl[0]=round(h*255.0);
    //     hsl[1]=round(s*255.0);
    //     hsl[2]=round(l*255.0);
    // };
    // inline void hsl_to_rgb(int rgb[3], int hsl[3]){
    //     double h=hsl[0]/255.0;
    //     double s=hsl[1]/255.0;
    //     double l=hsl[2]/255.0;
    //     int    f=int(h/60.0);
    //     double p=l*(1-s), q=l*(1-s*f);

    //     double r, g, b;
    //     int    rgb[3];
    //     if(f==0){
    //         r=q; g=p; b=p;
    //     }else if(f==1){
    //         r=p; g=q; b=p;
    //     }else if(f==2){
    //         r=p; g=q; b=p;
    //     }
    //     rgb[0]=round(r*255.0);
    //     rgb[1]=round(g*255.0);
    //     rgb[2]=round(b*255.0);
    // };
};

#endif