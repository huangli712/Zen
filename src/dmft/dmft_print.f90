!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : 
!!! source  : dmft_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           02/23/2021 by li huang (last modified)
!!! purpose : 
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dmft_print_header
!!
  subroutine dmft_print_header()
     use constants, only : mystd

     use version, only : V_FULL_IQ
     use version, only : V_AUTH_IQ
     use version, only : V_INST_IQ
     use version, only : V_MAIL_IQ
     use version, only : V_GPL3_IQ

     use control, only : cname
     use control, only : nprocs

     implicit none

! local variables
! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

# if defined (MPI)

     write(mystd,'(2X,a)') cname//' (Parallelized Edition)'

# else   /* MPI */

     write(mystd,'(2X,a)') cname//' (Sequential Edition)'

# endif  /* MPI */

     write(mystd,'(2X,a)') 'A Modern Continuous Time Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//V_FULL_IQ//' (built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//V_AUTH_IQ//' ('//V_INST_IQ//')'
     write(mystd,'(2X,a)') 'Support: '//V_MAIL_IQ
     write(mystd,'(2X,a)') 'License: '//V_GPL3_IQ
     write(mystd,*)

     write(mystd,'(2X,a)') 'start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', 1

# endif  /* MPI */

     return
  end subroutine dmft_print_header

!!
!! @sub: dmft_print_footer
!!
  subroutine dmft_print_footer()
     use constants, only : dp
     use constants, only : mystd

     use control, only : cname

     implicit none

! local variables
! string for current date and time
     character (len = 20) :: date_time_string

! used to record the time usage information
     real(dp) :: tot_time

! obtain time usage information
     call cpu_time(tot_time)

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') cname//' >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') cname//' >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') cname//' >>> happy ending at '//date_time_string

     return
  end subroutine dmft_print_footer
