#ifndef _VECTOR_MATH_H_
#define _VECTOR_MATH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "json_parser.h"

// Vector type
typedef double V3[3];

static inline double sqr(double v){
    return v*v;
}

static inline void v3_zero(V3 vector){
    vector[0] = 0;
    vector[1] = 0;
    vector[2] = 0;
}

static inline void v3_copy(V3 in, V3 out){
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
}

static inline void normalize(double *v) {
    double len = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    len = sqrt(len);
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
}


static inline double v3_len(V3 a){
    return sqrt(sqr(a[0])+sqr(a[1])+sqr(a[2]));
}


static inline void v3_add(V3 a, V3 b, V3 c){
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}


static inline void v3_subtract(V3 a, V3 b, V3 c){
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}


static inline void v3_scale(V3 a, double s, V3 b){
    b[0] = s * a[0];
    b[1] = s * a[1];
    b[2] = s * a[2];
}


static inline double v3_dot(V3 a, V3 b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


static inline void v3_cross(V3 a, V3 b, V3 c){
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}



void v3_reflect(V3 v, V3 n, V3 v_r){
    normalize(n);
    double scalar = 2.0 * v3_dot(n, v);
    V3 tmp_vector;
    v3_scale(n, scalar, tmp_vector);
    v3_subtract(v, tmp_vector, v_r);
}

#endif
