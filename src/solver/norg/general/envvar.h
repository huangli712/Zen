/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
    Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2021-2022
*/

#ifndef _ENVVAR_H_
#define _ENVVAR_H_

#ifdef _MSC_VER
#define _GHOST_HUNTING_
#else
#define _ASSERTIONFORVSCODE_
#endif

#ifdef _GHOST_HUNTING_
#define _CHECK_BOUNDS_
#define _CHECK_DIMENSION_MATCH_
#define _ASSERTION_
const int debug = 1;
#else
const int debug = 0;
#endif

#endif /* _ENVVAR_H_ */
