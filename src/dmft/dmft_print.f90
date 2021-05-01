!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_print_header
!!!           dmft_print_footer
!!!           dmft_print_summary
!!!           dmft_print_system
!!! source  : dmft_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           05/01/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dmft_print_header
!!
!! print the startup information for the dmft/jacaranda code
!!
  subroutine dmft_print_header()
     use constants, only : mystd

     use version, only : V_FULL_ZD
     use version, only : V_AUTH_ZD
     use version, only : V_INST_ZD
     use version, only : V_MAIL_ZD
     use version, only : V_GPL3_ZD

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

     write(mystd,'(2X,a)') 'A Modern Dynamical Mean-Field Theory Booster'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//V_FULL_ZD//' (built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//V_AUTH_ZD//' ('//V_INST_ZD//')'
     write(mystd,'(2X,a)') 'Support: '//V_MAIL_ZD
     write(mystd,'(2X,a)') 'License: '//V_GPL3_ZD
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
!! print the ending information for the dmft/jacaranda
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

!!
!! @sub dmft_print_summary
!!
!! print the running parameters, only for reference
!!
  subroutine dmft_print_summary()
     use constants, only : mystd

     use control ! ALL
     use context, only : qdim, qbnd

     implicit none

     write(mystd,'(2X,a)') '[configuration parameters] -> general'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,a10,  2X,a8)') 'model  / value :', trim(model) , 'type : i'

     write(mystd,'(2X,a)') '[configuration parameters] -> core control'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'task   / value :', task  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'axis   / value :', axis  , 'type : i'
     write(mystd,'(4X,a16,l10,  2X,a8)') 'lfermi / value :', lfermi, 'type : i'
     write(mystd,'(4X,a16,l10,  2X,a8)') 'ltetra / value :', ltetra, 'type : i'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'beta   / value :', beta  , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'mc     / value :', mc    , 'type : d'

     write(mystd,'(2X,a)') '[configuration parameters] -> kohn-sham data'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nsort  / value :', nsort , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'natom  / value :', natom , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nband  / value :', nband , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nkpt   / value :', nkpt  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nspin  / value :', nspin , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'ntet   / value :', ntet  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'ngrp   / value :', ngrp  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nwnd   / value :', nwnd  , 'type : i'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'scale  / value :', scale , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'fermi  / value :', fermi , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'volt   / value :', volt  , 'type : d'

     write(mystd,'(2X,a)') '[configuration parameters] -> local orbital projector'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'qdim   / value :', qdim  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'qbnd   / value :', qbnd  , 'type : i'

     write(mystd,'(2X,a)') '[configuration parameters] -> self-energy functions'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nsite  / value :', nsite , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nmesh  / value :', nmesh , 'type : i'

     write(mystd,*)

     return
  end subroutine dmft_print_summary

!!
!! @sub dmft_print_system
!!
!! print the system information, only for reference
!!
  subroutine dmft_print_system()
     use constants, only : mystd

     use control, only : nsort, natom 
     use control, only : nspin
     use control, only : ngrp, nwnd
     use control, only : nsite, nmesh

     use context, only : i_grp, i_wnd, g_imp, w_imp
     use context, only : shell, corr, site, l, ndim
     use context, only : bmin, bmax, nbnd
     use context, only : sorts, atoms, sortn, coord
     use context, only : sigdc, sig_l

     implicit none

! local variables
! loop index
     integer :: i
     integer :: s
     integer :: p

! print lattice structure
     write(mystd,'(2X,a)') "[system information] -> lattice    -> sorts" 
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,nsort
         write(mystd,'(4X,a6,i3,2X,a8,a3,i3)') "sort :", s, "symbol :", sorts(s), sortn(s)
     enddo ! over s={1,nsort} loop

     write(mystd,'(2X,a)') "[system information] -> lattice    -> atoms" 
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,natom
         write(mystd,'(4X,a6,i3,2X,a8,a3,1X,3f9.6)') "atom :", s, "symbol :", atoms(s), coord(s,:)
     enddo ! over s={1,natom} loop
     write(mystd,*)

! print quantum impurities
     write(mystd,'(2X,a)') "[system information] -> impurities -> sig_l" 
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,nsite
         do p=1,nspin
             write(mystd,'(4X,a10,i3,2X,a6,i3)') "impurity :", s, "spin :", p
             do i=1,ndim(i_grp(s))
                 write(mystd,'(4X,a1,i3,2f10.5)') ">", i, real(sig_l(1,i,i,p,s)), real(sig_l(nmesh,i,i,p,s))
             enddo ! over i={1,ndim(i_grp(s))} loop
         enddo ! over p={1,nspin} loop
     enddo ! over s={1,nsite} loop

     write(mystd,'(2X,a)') "[system information] -> impurities -> sigdc" 
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,nsite
         do p=1,nspin
             write(mystd,'(4X,a10,i3,2X,a6,i3)') "impurity :", s, "spin :", p
             do i=1,ndim(i_grp(s))
                 write(mystd,'(4X,a1,i3,f10.5)') ">", i, real(sigdc(i,i,p,s))
             enddo ! over i={1,ndim(i_grp(s))} loop
         enddo ! over p={1,nspin} loop
     enddo ! over s={1,nsite} loop

     write(mystd,'(2X,a)') "[system information] -> impurities -> mappings" 
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,nsite
         write(mystd,'(4X,a10,i3,2X,a8,i3)') "impurity :", s, "group  :", i_grp(s)
         write(mystd,'(4X,a10,i3,2X,a8,i3)') "impurity :", s, "window :", i_wnd(s)
     enddo ! over s={1,nsite} loop
     write(mystd,*)

! print local projectors
     write(mystd,'(2X,a)') "[system information] -> projectors -> groups"
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,ngrp
         write(mystd,'(4X,a8,i3)') "group  :", s
         write(mystd,'(4X,a9,2X,a)') "> shell :", shell(s)
         write(mystd,'(4X,a9,l3)') ">  corr :", corr(s)
         write(mystd,'(4X,a9,i3)') ">  site :", site(s)
         write(mystd,'(4X,a9,i3)') ">     l :", l(s)
         write(mystd,'(4X,a9,i3)') ">  ndim :", ndim(s)
     enddo ! over s={1,ngrp} loop

     write(mystd,'(2X,a)') "[system information] -> projectors -> windows"
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,nwnd
         write(mystd,'(4X,a8,i3)') "window :", s
         write(mystd,'(4X,a9,i3)') ">  bmin :", bmin(s)
         write(mystd,'(4X,a9,i3)') ">  bmax :", bmax(s)
         write(mystd,'(4X,a9,i3)') ">  nbnd :", nbnd(s)
     enddo ! over s={1,nwnd} loop

     write(mystd,'(2X,a)') "[system information] -> projectors -> mappings"
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     do s=1,ngrp
         write(mystd,'(4X,a8,i3,2X,a10,i3)') "group  :", s, "impurity :", g_imp(s)
     enddo ! over s={1,ngrp} loop
     do s=1,nwnd
         write(mystd,'(4X,a8,i3,2X,a10,i3)') "window :", s, "impurity :", w_imp(s)
     enddo ! over s={1,nwnd} loop

     write(mystd,*)

     return
  end subroutine dmft_print_system
