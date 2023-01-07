/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _MYMKL_H_
#define _MYMKL_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

#define _LAPACKE_

#include "miblas.h"

#ifdef _LAPACKE_
#include "milapacke.h"
#else
#include "milapack.h"
#endif


#endif /* _MYMKL_H_ */
