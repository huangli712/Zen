This module wraps the most useful message passing interface (MPI) calls by using generic programming techniques. It supports most of the collective operations (such as BCAST, GATHER, REDUCE, etc.). These subroutines are

* MPI\_INIT()
* MPI\_FINALIZE()
* MPI\_WTIME()
* MPI\_WTICK()
* MPI\_BARRIER()
* MPI\_DIMS\_CREATE()
* MPI\_CART\_CREATE()
* MPI\_CART\_COORDS()
* MPI\_COMM\_SPLIT()
* MPI\_COMM\_RANK()
* MPI\_COMM\_SIZE()
* MPI\_GET\_PROCESSOR_NAME()
* MPI\_BCAST()
* MPI\_GATHER()
* MPI\_GATHERV()
* MPI\_ALLGATHER()
* MPI\_ALLGATHERV()
* MPI\_REDUCE()
* MPI\_ALLREDUCE()

etc. However, none of the point-to-point operations is supported. In the module, we also try to implement a light-weight error handler. Enjoy it!

## Type

module

## Source

`src/m_mpi.f90`

## Usage

**(1)** Include mpi support.

```fortran
use mmpi
```

Please pay attention to the module name. It is `mmpi`, instead of `mpi`.

**(2)** Init mpi environment.

```fortran
call mp_init()            ! init mpi environment
call mp_comm_rank(myid)   ! get current process it
call mp_comm_size(nprocs) ! get number of processes
```

**(3)** Broadcast data.

```fortran
real(dp), allocatable :: real_data(:,:,:)
integer, allocatable :: int_data(:)
complex(dp), allocatable :: cmplx_data(:,:,:,:)

call mp_bcast(real_data, master)
call mp_bcast(int_data, master)
call mp_bcast(cmplx_data, master)
```

Here master == 0 which means the master node/root process.

**(4)** Allreduce data.

```fortran
real(dp), allocatable :: real_data(:)
real(dp), allocatable :: real_data_mpi(:)

call mp_allreduce(real_data, real_data_mpi)     ! all-readuce data
real_data = real_data_mpi / number_of_processes ! calculate the average
```

**(5)** Setup barrier.

```fortran
call mp_barrier()
```

**(6)** Finialize mpi environment.

```fortran
call mp_finalize()
```
