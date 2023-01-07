/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2012 - 2017
*/

#ifndef _EIGSYM_H_
#define _EIGSYM_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

Int eigsym(MatReal &z, VecReal &d, Int ifv = 1);
void tred2(MatReal &z, VecReal &d, VecReal &e, Int ifv = 1);
Int tqli(MatReal &z, VecReal &d, VecReal &e, Int ifv = 1);
Int tqlis(VecReal &z, VecReal &d, VecReal &e, Int ifv = 1);

#endif /* _EIGSYM_H_ */
