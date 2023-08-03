#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define CubiCCC_IMPORT
#include "CubiCCC.h"

static double * init_coeff_matrix(double *, uint32_t, uint32_t);
static double f_kernel(double);
static double * f_cubic_convolution(double *, uint32_t, uint32_t, double *, uint32_t, uint32_t);

double * CubiCCC_interpolate(double * input_matrix, uint32_t input_width, uint32_t input_height, double * output_matrix, uint32_t output_width, uint32_t output_height){
    if((output_width == 0) || (output_height == 0)) return NULL;

    //if((input_width < 3) || (input_height < 3)) do something else
    double * coeff_matrix = init_coeff_matrix(input_matrix, input_width, input_height);

    //Input and output imagined as rectangles from (0, 0) to (width - 1, height - 1) hence width/height horizontal/vertical points counting (x, 0) and (0, y)
    //Spacing between nodes = 1 for input matrix

    //double in_to_out[2] = {(double) (input_width - 1)/(output_width - 1), (double) (input_height - 1)/(output_height - 1)};
    f_cubic_convolution(coeff_matrix, input_width + 2, input_height + 2, output_matrix, output_width, output_height);
    //free(output_matrix);
    free(coeff_matrix);
    return output_matrix;
    //return coeff_matrix;
}

static double * init_coeff_matrix(double * input, uint32_t width, uint32_t height){
    uint32_t i, j;
    //c_m = coefficient matrix
    uint32_t c_m_width = width + 2;
    uint32_t c_m_height = height + 2;
    double * c_m = (double *) calloc(c_m_width*c_m_height, sizeof(double));

    for(i = 1; i < c_m_height - 1; i++)
        for(j = 1; j < c_m_width - 1; j++)
            c_m[(i * c_m_width) + j] = input[((i - 1) * width) + (j - 1)];

    for(i = 1; i < c_m_height - 1; i++){
        if(width > 2){
            c_m[(i * c_m_width)] = (3 * input[((i - 1) * width)]) - (3 * input[((i - 1) * width) + 1]) + (input[((i - 1) * width) + 2]);
            c_m[(i * c_m_width) + (c_m_width - 1)] = (3 * input[((i - 1) * width) + (width - 1)]) - (3 * input[((i - 1) * width) + (width - 2)]) + (input[((i - 1) * width) + (width - 3)]);
        }
    }
    for(j = 1; j < c_m_width - 1; j++){
        if(height > 2){
            c_m[j] = (3 * input[(j - 1)]) - (3 * input[width + (j - 1)]) + (input[(2 * width) + (j - 1)]);
            c_m[((c_m_height - 1) * c_m_width) + j] = (3 * input[((height - 1) * width) + (j - 1)]) - (3 * input[((height - 2) * width) + (j - 1)]) + (input[((height - 3) * width) + (j - 1)]);
        }
    }
    if(width > 2){
        c_m[0] = (3 * c_m[1]) - (3 * c_m[2]) + c_m[3];
        c_m[(c_m_width - 1)] = (3 * c_m[(c_m_width - 2)]) - (3 * c_m[(c_m_width - 3)]) + c_m[(c_m_width - 4)];
        c_m[((c_m_height - 1) * c_m_width)] = (3 * c_m[((c_m_height - 1) * c_m_width) + 1]) - (3 * c_m[((c_m_height - 1) * c_m_width) + 2]) + c_m[((c_m_height - 1) * c_m_width) + 3];
        c_m[((c_m_height - 1) * c_m_width) + (c_m_width - 1)] = (3 * c_m[((c_m_height - 1) * c_m_width) + (c_m_width - 2)]) - (3 * c_m[((c_m_height - 1) * c_m_width) + (c_m_width - 3)]) + c_m[((c_m_height - 1) * c_m_width) + (c_m_width - 4)];
    }
    return c_m;
}

static double f_kernel(double s){
    if(s >= 0){
        if(s < 1) return (1.5 * s*s*s) - (2.5 * s*s) + 1;
        if(s < 2) return 2.5*(s*s) + 2 - 0.5*(s*s*s) - 4*s;
        return 0;
    } else {
        return f_kernel(-s);
    }
}

double * f_cubic_convolution(double * coeff_matrix, uint32_t c_m_width, uint32_t c_m_height, double * output_matrix, uint32_t output_width, uint32_t output_height){
    uint32_t i, j;
    int8_t k, l;
    double ratio_x = (double) (c_m_width - 3)/(output_width - 1);
    double ratio_y = (double) (c_m_height - 3)/(output_height - 1);
    double * kernel_x = (double *) calloc(4*output_width, sizeof(double));
    double * kernel_y = (double *) calloc(4*output_height, sizeof(double));
    double x, y;
    double floor_x, floor_y;

    if(output_width > 1){
        for(i = 0; i < output_width; i++){
            x = (double) i * ratio_x;
            floor_x = floor(x);
            for(k = -1; k < 3; k++)
                kernel_x[(i * 4) + (k + 1)] = f_kernel(x - (floor_x + k));
        }
    } else {
        kernel_x[1] = 1;
    }
    if(output_height > 1){
        for(i = 0; i < output_height; i++){
            y = (double) i * ratio_y;
            floor_y = floor(y);
            for(k = -1; k < 3; k++)
                kernel_y[(i * 4) + (k + 1)] = f_kernel(y - (floor_y + k));
        }
    } else {
        kernel_y[1] = 1;
    }

    for(i = 0; i < output_height; i++){
        for(j = 0; j < output_width; j++){
            output_matrix[(i * output_width) + j] = 0;
            y = (double) i * ratio_y;
            floor_y = floor(y);
            x = (double) j * ratio_x;
            floor_x = floor(x);
            for(k = -1; k < 3; k++)
                for(l = -1; l < 3; l++)
                    output_matrix[(i * output_width) + j] += coeff_matrix[(((uint32_t)floor_y + 1 + l) * c_m_width) + ((uint32_t)floor_x + 1 + k)] * kernel_x[(j * 4) + (k + 1)] * kernel_y[(i * 4) + (l + 1)];
        }
    }

    free(kernel_x);
    free(kernel_y);

    return output_matrix;
}
