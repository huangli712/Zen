#pragma once

/*
Linear equation solution by Gauss-Jordan elimination

find matrix x so that a x = b
where a is a square matrix
and b is a rectangular matrix

Tranlated from Numerical Recipies in c++ (Third Edition)
by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China), 2021-02-18
*/

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

// find matrix x so that a x = b, 
// where a is a square matrix, 
// and b is a rectangular matrix
// on output, a is replaced by its matrix inverse, 
// and b is replaced by the corresponding set of solution vectors
// return 1 if a is singular
// return 0 if meaningful results are returned in a and b
Int gaussj(MatReal& a, MatReal& b);

// on output, a is replaced by its matrix inverse
// return 1 if a is singular
// return 0 if a is inverted successfully
Int gaussj(MatReal& a);
