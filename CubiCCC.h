#ifndef CubiCCC_H
#define CubiCCC_H

#include <stdint.h>

#ifdef CubiCCC_IMPORT
    #define EXTERN
#else
    #define EXTERN extern
#endif // CubiCCC_IMPORT

//Arguments:    input matrix, ncols, nrows, output matrix, ncols, nrows
//Returns:      output matrix
EXTERN double * CubiCCC_interpolate(double *, uint32_t, uint32_t, double *, uint32_t, uint32_t);


#undef CubiCCC_IMPORT
#undef EXTERN
#endif // CubiCCC_H
