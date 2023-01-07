#pragma once

/*
code developed by
    Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/
#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

/*---------------------------------------------------------------------------   */
/*  Example program for solving symmetric positive definite system of equations.*/
/*  Medium case: make use of the user-defined stopping test.                    */
/*---------------------------------------------------------------------------   */

// NOT Tested yet!! --2022.7.19
template<typename F>
int conjugate_gradient(F& A, VecReal& x, const VecReal& b, MKL_INT expected_itercount = 100, Int max_itercount = 10000, Int test_mode = 0)
    //	Conjugate Gradient method to solve the x in form of A·x = b, which A is a
    //  symmetric positive definite matrix. In here the A is suggest as the sparse
    //  matrix form.
    //Input:
    //  x           is the initial guess vector, first we input the expected soluction.
    //  A           typename F(matrix) A needed with well define operator*(which to returns A·x)
    //  b           is the aim of vector of A·x = b.
    
{
#ifdef _CHECK_DIMENSION_MATCH_
    ASSERT_EQ(x.size(), b.size());
#endif
    const MKL_INT n = b.size();
    MKL_INT rci_request, itercount, i;
    double rhs[n];
    VecReal rhs_v(n, &rhs[0]);
    /*---------------------------------------------------------------------------*/
    /* Allocate storage for the solver ?par and temporary storage tmp            */
    /*---------------------------------------------------------------------------*/
    MKL_INT length = 128;
    MKL_INT ipar[128];
    double dpar[128], tmp[4 * n];
    /*---------------------------------------------------------------------------*/
    /* Some additional variables to use with the RCI (P)CG solver                */
    /*---------------------------------------------------------------------------*/
    double solution[n];
    VecReal expected_sol(b);
    VecReal x_result(n, &solution[0]);
    double euclidean_norm, temp[n];
    VecReal temp_v(n, &temp[0]);
    double eone = -1.E0;
    MKL_INT ione = 1;

    //DBG("To here is right!" );

    /*---------------------------------------------------------------------------*/
    /* Initialize the right hand side through matrix-vector product              */
    /*---------------------------------------------------------------------------*/
    //mkl_dcsrsymv(&tr, &n, a, ia, ja, expected_sol, rhs);
    rhs_v = A * x;
    //DBG("To here is right!");
    /*---------------------------------------------------------------------------*/
    /* Initialize the initial guess                                              */
    /*---------------------------------------------------------------------------*/
    for (i = 0; i < n; i++)
    {
        solution[i] = 0.E0;
    }

    //DBG("To here is right!");
    /*---------------------------------------------------------------------------*/
    /* Initialize the solver                                                     */
    /*---------------------------------------------------------------------------*/
    dcg_init(&n, solution, rhs, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0)
        goto failure;
    /*---------------------------------------------------------------------------*/
    /* Set the desired parameters:                                               */
    /* INTEGER parameters:                                                        */
    /* set the maximal number of iterations to 100                               */
    /* LOGICAL parameters:                                                       */
    /* run the Preconditioned version of RCI (P)CG with preconditioner C_inverse */
    /* DOUBLE parameters                                                         */
    /* -                                                                         */
    /*---------------------------------------------------------------------------*/
    ipar[4] = max_itercount;
    //DBG("To here is right!");
    /*---------------------------------------------------------------------------*/
    /* Check the correctness and consistency of the newly set parameters         */
    /*---------------------------------------------------------------------------*/
    dcg_check(&n, solution, rhs, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0)
        goto failure;
    /*---------------------------------------------------------------------------*/
    /* Compute the solution by RCI (P)CG solver                                  */
    /* Reverse Communications starts here                                        */
    /*---------------------------------------------------------------------------*/
rci:dcg(&n, solution, rhs, &rci_request, ipar, dpar, tmp);
    /*---------------------------------------------------------------------------*/
    /* If rci_request=0, then the solution was found according to the requested  */
    /* stopping tests. In this case, this means that it was found after 100      */
    /* iterations.                                                               */
    /*---------------------------------------------------------------------------*/
    if (rci_request == 0)
        goto getsln;
    /*---------------------------------------------------------------------------*/
    /* If rci_request=1, then compute the vector A*tmp[0]                      */
    /* and put the result in vector tmp[n]                                     */
    /*---------------------------------------------------------------------------*/
    if (rci_request == 1)
    {
        //mkl_dcsrsymv(&tr, &n, a, ia, ja, tmp, &tmp[n]);//
        VecReal tmp_i(n, &tmp[0]), tmp_ii(n, &tmp[n]);
        tmp_ii = A * x;
        goto rci;
    }
    /*---------------------------------------------------------------------------*/
    /* If rci_request=2, then do the user-defined stopping test: compute the     */
    /* Euclidean norm of the actual residual using MKL routines and check if     */
    /* it is less than 1.E-8                                                     */
    /*---------------------------------------------------------------------------*/
    if (rci_request == 2)
    {
        //mkl_dcsrsymv(&tr, &n, a, ia, ja, solution, temp);
        temp_v = A * x;
        daxpy(&n, &eone, rhs, &ione, temp, &ione);
        euclidean_norm = dnrm2(&n, temp, &ione);
        /*---------------------------------------------------------------------------*/
        /* The solution has not been found yet according to the user-defined stopping */
        /* test. Continue RCI (P)CG iterations.                                      */
        /*---------------------------------------------------------------------------*/
        if (euclidean_norm > 1.e-8)
            goto rci;
        /*---------------------------------------------------------------------------*/
        /* The solution has been found according to the user-defined stopping test   */
        /*---------------------------------------------------------------------------*/
        else
            goto getsln;
    }
    /*---------------------------------------------------------------------------*/
    /* If rci_request=anything else, then dcg subroutine failed                  */
    /* to compute the solution vector: solution[n]                               */
    /*---------------------------------------------------------------------------*/
    goto failure;
    /*---------------------------------------------------------------------------*/
    /* Reverse Communication ends here                                           */
    /* Get the current iteration number into itercount                           */
    /*---------------------------------------------------------------------------*/
getsln:dcg_get(&n, solution, rhs, &rci_request, ipar, dpar, tmp, &itercount);
    /*---------------------------------------------------------------------------*/
    /* Print solution vector: solution[n] and number of iterations: itercount    */
    /*---------------------------------------------------------------------------*/
    x = x_result;
    if(test_mode){
        printf("The system has been solved\n");
        printf("The following solution obtained\n");
        for (i = 0; i < n / 2; i++)
            printf("%6.3f  ", solution[i]);
        printf("\n");
        for (i = n / 2; i < n; i++)
            printf("%6.3f  ", solution[i]);
        printf("\nExpected solution is\n");
        for (i = 0; i < n / 2; i++)
        {
            printf("%6.3f  ", expected_sol[i]);
            expected_sol[i] -= solution[i];
        }
        printf("\n");
        for (i = n / 2; i < n; i++)
        {
            printf("%6.3f  ", expected_sol[i]);
            expected_sol[i] -= solution[i];
        }
        printf("\nNumber of iterations: %d\n", itercount);
    }
    i = 1;
    euclidean_norm = dnrm2(&n, expected_sol.begin(), &i);

    /*-------------------------------------------------------------------------*/
    /* Release internal MKL memory that might be used for computations         */
    /* NOTE: It is important to call the routine below to avoid memory leaks   */
    /* unless you disable MKL Memory Manager                                   */
    /*-------------------------------------------------------------------------*/
    MKL_Free_Buffers();

    if (itercount == expected_itercount && euclidean_norm < 1.0e-12)
    {
        printf("This example has successfully PASSED through all steps of computation!\n");
        return 0;
    }
    else
    {
        printf("This example may have FAILED as either the number of iterations differs\n");
        printf("from the expected number of iterations %d, or the ", expected_itercount);
        printf("computed solution\ndiffers much from the expected solution ");
        printf("(Euclidean norm is %e), or both.\n", euclidean_norm);
        return 1;
    }
    /*-------------------------------------------------------------------------*/
    /* Release internal MKL memory that might be used for computations         */
    /* NOTE: It is important to call the routine below to avoid memory leaks   */
    /* unless you disable MKL Memory Manager                                   */
    /*-------------------------------------------------------------------------*/
failure:printf("This example FAILED as the solver has returned the ERROR code %d", rci_request);
    MKL_Free_Buffers();
    return 1;
  
}

template<typename F>
int conjugate_gradient_simple(F& A, VecReal& x, const VecReal& b, MKL_INT expected_itercount = 100, Int max_itercount = 10000, Int test_mode = 0)
{
    VecReal r(b - A * x);
    VecReal p(r);
    Real rsold(DOT(r, r));
    VecReal ap(A * p);
    Real alpha(rsold / DOT(p, ap));
    x = x + alpha * p;
    r = r - alpha * ap;
    Real rsnew(DOT(r, r));

    for_Int(i, 0, b.size()) {
        if (sqrt(rsnew) < 1e-10)break;
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
        ap = A * p;
        alpha = rsold / DOT(p, ap);
        x = x + alpha * p;
        r = r - alpha * ap;
        rsnew = DOT(r, r);
    }
    return 0;

//function x = conjgrad(A, b, x)
//    r = b - A * x;
//    p = r;
//    rsold = r' * r;
//
//    for i = 1:length(b)
//        Ap = A * p;
//        alpha = rsold / (p' * Ap);
//        x = x + alpha * p;
//        r = r - alpha * Ap;
//        rsnew = r' * r;
//        if sqrt(rsnew) < 1e-10
//              break
//        end
//        p = r + (rsnew / rsold) * p;
//        rsold = rsnew;
//    end
//end

}