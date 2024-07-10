#include "GRAPH.h"

void image_pixels(const char* png_path, unsigned char *pixels, int numpx, int numpy)
{	
    png_structp png_ptr;  
    png_infop info_ptr;  
    FILE *fp=fopen(png_path, "wb");  
    if(!fp){
        printf("[ERROR] Unable to create %s using fopen.\n", png_path);
        exit(EXIT_FAILURE);
    }
    png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);  
    if(png_ptr==NULL){  
        printf("[ERROR] Unable to create %s using png_create_write_struct.\n", png_path);
        fclose(fp);
        exit(EXIT_FAILURE);
    }  
    info_ptr=png_create_info_struct(png_ptr);
    if(info_ptr==NULL){  
        printf("[ERROR] Unable to create %s using png_create_info_struct.\n", png_path);
        png_destroy_write_struct(&png_ptr, NULL);  
        exit(EXIT_FAILURE);
    }  
    png_init_io(png_ptr, fp);  
    png_set_IHDR(png_ptr, info_ptr, numpx, numpy, 8, 
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE); 
    png_colorp palette=(png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH*sizeof(png_color));
    if(!palette){
        printf("[ERROR] Unable to create palette using png_malloc.\n");
        fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        exit(EXIT_FAILURE);
    }
    png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);  
    png_write_info(png_ptr, info_ptr);  
    png_set_packing(png_ptr);
    png_bytepp rows=(png_bytepp)png_malloc(png_ptr, numpy*sizeof(png_bytep));
    for(int i=0;i<numpy;++i){
        rows[numpy-i-1]=(png_bytep)(pixels+(i)*numpx*3);
    }
    png_write_image(png_ptr, rows);  
    delete[] rows;
    png_write_end(png_ptr, info_ptr);  
    png_free(png_ptr, palette);  
    palette=NULL;  
    png_destroy_write_struct(&png_ptr, &info_ptr);  
    fclose(fp);   
}

GRAPH::GRAPH(double width, double height, int resolution)
{
    image.pwidth=round(width*resolution); 
    image.pheight=round(height*resolution);
    callocate_2d(&image.pixels, image.pheight, image.pwidth, image.black);
}

GRAPH::~GRAPH()
{
    deallocate_2d(image.pixels, image.pheight);
    deallocate(points);
}

void GRAPH::set_xlim(double xmin, double xmax)
{
    area.start_x=xmin; 
    area.end_x=xmax;
    area.spacing_x=xmax-xmin;
}

void GRAPH::set_ylim(double ymin, double ymax)
{
    area.start_y=ymin; 
    area.end_y=ymax;
    area.spacing_y=ymax-ymin;
}

void GRAPH::set_xtick_spacing(double major_spacing, double minor_spacing)
{
    top.major_tick_spacing=major_spacing;
    down.major_tick_spacing=major_spacing;
    top.minor_tick_spacing=minor_spacing;
    down.minor_tick_spacing=minor_spacing;
}

void GRAPH::set_ytick_spacing(double major_spacing, double minor_spacing)
{
    left.major_tick_spacing=major_spacing;
    right.major_tick_spacing=major_spacing;
    left.minor_tick_spacing=minor_spacing;
    right.minor_tick_spacing=minor_spacing;
}

void GRAPH::set_tick_in(bool is_tick_in)
{
    left.is_tick_in=is_tick_in;
    right.is_tick_in=is_tick_in;
    top.is_tick_in=is_tick_in;
    down.is_tick_in=is_tick_in;
}

void GRAPH::set_top_visible(bool is_visible)
{
    top.is_visible=is_visible;
}

void GRAPH::set_right_visible(bool is_visible)
{
    right.is_visible=is_visible;
}

void GRAPH::draw_xaxis(AXIS *axis, double x1, double x2, double y)
{
    if(!axis->is_visible) return;
    POINT p1, p2;
    p1.black=p2.black=axis->black;
    p1.style=p2.style='s';
    p1.psize=p2.psize=axis->line_pwidth;
    p1.is_in_area_checked=p2.is_in_area_checked=axis->is_in_area_checked;
    p1.x=x1; p2.x=x2; p1.y=p2.y=y;
    draw_line(&p1, &p2);
    double x_major_tick_length=axis->major_tick_psize*area.ratio_y;
    double x_minor_tick_length=axis->minor_tick_psize*area.ratio_y;
    int    n_major_tick=round(area.spacing_x/axis->major_tick_spacing);
    int    n_minor_tick=round(axis->major_tick_spacing/axis->minor_tick_spacing);
    int    dr=axis->is_tick_in?1:-1;
    for(int i=0;i<n_major_tick;i++){
        p1.x=x1+i*axis->major_tick_spacing;
        p2.x=p1.x; p2.y=y+x_major_tick_length*dr;
        draw_line(&p1, &p2);
        for(int j=0;j<n_minor_tick;j++){
            p1.x=p1.x+j*axis->minor_tick_spacing;
            p2.x=p1.x; p2.y=y+x_minor_tick_length*dr;
            draw_line(&p1, &p2);
        }
    }
    p1.x=x2;
    p2.x=p1.x; p2.y=y+x_major_tick_length*dr;
    draw_line(&p1, &p2);
}

void GRAPH::draw_yaxis(AXIS *axis, double y1, double y2, double x)
{
    if(!axis->is_visible) return;
    POINT p1, p2;
    p1.black=p2.black=axis->black;
    p1.style=p2.style='s';
    p1.psize=p2.psize=axis->line_pwidth;
    p1.is_in_area_checked=p2.is_in_area_checked=axis->is_in_area_checked;
    p1.x=x; p1.y=y1;
    p2.x=x; p2.y=y2;
    draw_line(&p1, &p2);
    double y_major_tick_length=axis->major_tick_psize*area.ratio_x;
    double y_minor_tick_length=axis->minor_tick_psize*area.ratio_x;
    int    n_major_tick=round(area.spacing_y/axis->major_tick_spacing);
    int    n_minor_tick=round(axis->major_tick_spacing/axis->minor_tick_spacing);
    int    dr=axis->is_tick_in?1:-1;
    for(int i=0;i<n_major_tick;i++){
        p1.y=y1+i*axis->major_tick_spacing;
        p2.x=x+y_major_tick_length*dr; p2.y=p1.y;
        draw_line(&p1, &p2); 
        for(int j=0;j<n_minor_tick;j++){
            p1.y=p1.y+j*axis->minor_tick_spacing;
            p2.x=x+y_minor_tick_length*dr; p2.y=p1.y;
            draw_line(&p1, &p2);
        }
    }
    p1.y=y2;
    p2.x=x+y_major_tick_length*dr; p2.y=p1.y;
    draw_line(&p1, &p2); 
}

void GRAPH::scatter(double *x, double *y, double *value, int num)
{
    strcpy(style, "scatter");
    nump=num;
    mallocate(&points, num);
    double max_value=0.0;
    for(int i=0;i<num;i++){
        points[i].x=x[i];
        points[i].y=y[i];
        if(max_value<value[i]){
            max_value=value[i];
        }
    }
    for(int i=0;i<num;i++){
        points[i].black=value[i]/max_value;
    }
}

void GRAPH::line(double *x, double *y, int num)
{
    strcpy(style, "line");
    nump=num;
    mallocate(&points, num);
    for(int i=0;i<num;i++){
        points[i].x=x[i];
        points[i].y=y[i];
    }
}

void GRAPH::draw(const char *png_path)
{
    area.start_px=round(area.border_width_factor*image.pwidth);
    area.start_py=round(area.border_height_factor*image.pheight);
    area.pwidth=round((1.0-area.border_width_factor*2.0)*image.pwidth); 
    area.pheight=round((1.0-area.border_height_factor*2.0)*image.pheight);
    area.end_px=area.start_px+area.pwidth;
    area.end_py=area.start_py+area.pheight;
    area.ratio_x=area.spacing_x/area.pwidth;
    area.ratio_y=area.spacing_y/area.pheight;
    if(0==strcmp(style, "scatter")){
        for(int i=0;i<nump;i++){
            draw_marker(&points[i]);
        }
    }else if(0==strcmp(style, "line")){
        for(int i=1;i<nump;i++){
            draw_line(&points[i-1], &points[i]);
        }
    }else{
        printf("[ERROR] Unable to draw graph without data input");
        exit(EXIT_FAILURE);
    }
    draw_xaxis(&top, area.start_x, area.end_x, area.end_y);
    draw_xaxis(&down, area.start_x, area.end_x, area.start_y);
    draw_yaxis(&left, area.start_y, area.end_y, area.start_x);
    draw_yaxis(&right, area.start_y, area.end_y, area.end_x);
    int npixel=image.pwidth*image.pheight;
    double *wpixels; callocate(&wpixels, npixel, 0.0);
    for(int i=0;i<image.pheight;i++){
        for(int j=0;j<image.pwidth;j++){
            wpixels[i*image.pwidth+j]=image.pixels[i][j];
        }
    }
    unsigned char *pixels; mallocate(&pixels, 3*npixel);
    for(int i=0;i<npixel;i++){
        pixels[i*3]=pixels[i*3+1]=pixels[i*3+2]=round((1.0-wpixels[i])*255.0);
    }
    image_pixels(png_path, pixels, image.pwidth, image.pheight);
}

void GRAPH::get_pixel(int &i, int &j, double x, double y)
{
    i=area.start_px+round((x-area.start_x)/area.ratio_x);
    j=area.start_py+round((y-area.start_y)/area.ratio_y);
}

void GRAPH::draw_marker(POINT *point)
{
    double sizex=point->psize*area.ratio_x;
    double sizey=point->psize*area.ratio_y;
    double start_x=point->x-sizex/2.0, start_y=point->y-sizey/2.0;
    double end_x=point->x+sizex/2.0, end_y=point->y+sizey/2.0;
    int start_i, start_j, end_i, end_j;
    get_pixel(start_i, start_j, start_x, start_y);
    get_pixel(end_i, end_j, end_x, end_y);
    int dii=end_i-start_i, djj=end_j-start_j;
    int istep=dii>0?1:-1, jstep=djj>0?1:-1;
    switch(point->style)
    {
    case 's':
        for(int i=start_i-1;i<end_i;i+=istep){
            for(int j=start_j-1;j<end_j;j+=jstep){
                if(((!point->is_in_area_checked)&&i>=0&&i<image.pwidth&&j>=0&&j<image.pheight)||(i>=area.start_px-1&&i<area.end_px&&j>=area.start_py&&j<area.end_py)){
                    image.pixels[j][i]=point->black;
                }
            }
        }
        break;
    case 'c':
        {
            int num=0, sum_i=0, sum_j=0;
            for(int i=start_i-1;i<end_i;i+=istep){
                for(int j=start_j-1;j<end_j;j+=jstep){
                    num++;
                    sum_i+=i; sum_j+=j;
                }
            }
            double center_i=double(sum_i)/double(num);
            double center_j=double(sum_j)/double(num);
            double r2=point->psize*point->psize/4.0;
            for(int i=start_i-1;i<end_i;i+=istep){
                for(int j=start_j-1;j<end_j;j+=jstep){
                    if((i-center_i)*(i-center_i)+(j-center_j)*(j-center_j)<r2){
                        if(((!point->is_in_area_checked)&&i>=0&&i<image.pwidth&&j>=0&&j<image.pheight)||(i>=area.start_px-1&&i<area.end_px&&j>=area.start_py&&j<area.end_py)){
                            image.pixels[j][i]=point->black;
                        }
                    }
                }
            }
        }
        break;
    default:
        printf("[ERROR] Unrecognized marker style %c", point->style);
        exit(EXIT_FAILURE);
    }
}

void GRAPH::draw_line(POINT *start_point, POINT *end_point)
{
    if(fabs(start_point->x-end_point->x)<fabs(start_point->y-end_point->y)){
        double dy=area.ratio_y;
        POINT  p={start_point->x, start_point->y, start_point->black, start_point->style, start_point->psize, start_point->is_in_area_checked};
        double k=(start_point->x-end_point->x)/(start_point->y-end_point->y);
        double b=end_point->x-k*end_point->y;
        if(start_point->y<end_point->y){
            for(p.y=start_point->y;p.y<end_point->y;p.y+=dy){
                p.x=k*p.y+b;
                draw_marker(&p);
            }
        }else{
            for(p.y=start_point->y;p.y>end_point->y;p.y-=dy){
                p.x=k*p.y+b;
                draw_marker(&p);
            } 
        }
        draw_marker(end_point);
    }else{
        double dx=area.ratio_x;
        POINT  p={start_point->x, start_point->y, start_point->black, start_point->style, start_point->psize, start_point->is_in_area_checked};
        double k=(start_point->y-end_point->y)/(start_point->x-end_point->x);
        double b=end_point->y-k*end_point->x;
        if(start_point->x<end_point->x){
            for(p.x=start_point->x;p.x<end_point->x;p.x+=dx){
                p.y=k*p.x+b;
                draw_marker(&p);
            }
        }else{
            for(p.x=start_point->x;p.x>end_point->x;p.x-=dx){
                p.y=k*p.x+b;
                draw_marker(&p);
            } 
        }
        draw_marker(end_point);
    }
}