!!!-----------------------------------------------------------------------
!!! project : flink @ sakura
!!! program : version
!!! source  : m_version.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/26/2017 by li huang (created)
!!!           04/15/2019 by li huang (last modified)
!!! purpose : the purpose of this module is to define version strings.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! It is a common module which defines the current version of iQIST (and
!! the other numerical applications).
!!
!! List of suffixes:
!!    _IQ => iQIST
!!    _FL => FLINK
!!    _DF => DFAPP
!!
!! Usage
!! =====
!!
!! use version, only : V_FULL_IQ
!! implicit none
!!
!! print *, V_FULL_IQ
!!
!!

  module version
     implicit none

!!
!! @var V_FULL
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: V_FULL_IQ = 'v0.8.0 @ 2019.04.10D'
     character(len=20), public, parameter :: V_FULL_FL = 'v0.1.1 @ 2019.04.12D'
     character(len=20), public, parameter :: V_FULL_DF = 'v0.0.1 @ 2019.04.12D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR_IQ = 'v0.8.0'
     character(len=06), public, parameter :: V_CURR_FL = 'v0.1.1'
     character(len=06), public, parameter :: V_CURR_DF = 'v0.0.1'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE_IQ = '2019.04.10'
     character(len=11), public, parameter :: V_DATE_FL = '2019.04.12'
     character(len=11), public, parameter :: V_DATE_DF = '2019.04.12'

!!
!! @var V_STAT
!!
!! version string, only status info., D means devel, T testing, R released.
!!
     character(len=01), public, parameter :: V_STAT_IQ = 'D'
     character(len=01), public, parameter :: V_STAT_FL = 'D'
     character(len=01), public, parameter :: V_STAT_DF = 'D'

!!
!! @var V_AUTH
!!
!! version string, author info.
!!
     character(len=11), public, parameter :: V_AUTH_IQ = 'by li huang'
     character(len=11), public, parameter :: V_AUTH_FL = 'by li huang'
     character(len=11), public, parameter :: V_AUTH_DF = 'by li huang'

!!
!! @var V_INST
!!
!! version string, affiliation info.
!!
     character(len=36), public, parameter :: V_INST_IQ = 'China Academy of Engineering Physics'
     character(len=36), public, parameter :: V_INST_FL = 'China Academy of Engineering Physics'
     character(len=36), public, parameter :: V_INST_DF = 'China Academy of Engineering Physics'

!!
!! @var V_MAIL
!!
!! version string, email info.
!!
     character(len=22), public, parameter :: V_MAIL_IQ = 'lihuang.dmft@gmail.com'
     character(len=22), public, parameter :: V_MAIL_FL = 'lihuang.dmft@gmail.com'
     character(len=22), public, parameter :: V_MAIL_DF = 'lihuang.dmft@gmail.com'

!!
!! @var V_GPL3
!!
!! version string, license info.
!!
     character(len=36), public, parameter :: V_GPL3_IQ = 'GNU General Public License version 3'
     character(len=36), public, parameter :: V_GPL3_FL = 'GNU General Public License version 3'
     character(len=36), public, parameter :: V_GPL3_DF = 'GNU General Public License version 3'

  end module version
