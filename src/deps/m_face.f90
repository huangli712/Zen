!!!-----------------------------------------------------------------------
!!! project : flink @ sakura
!!! program : face
!!! source  : m_face.f90
!!! type    : module
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 12/31/2024 by li huang (created)
!!!           01/16/2025 by li huang (last modified)
!!! purpose : to support colorful outputs via ascii escape sequences.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module face
     use, intrinsic :: iso_fortran_env, only: int32

     implicit none

!!========================================================================
!!>>> declare local parameters                                         <<<
!!========================================================================

!! module parameters

     ! general parameters
     !
     ! upper case alphabet
     character(26), private, parameter :: UPPER_ALPHABET='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
     !
     ! lower case alphabet
     character(26), private, parameter :: LOWER_ALPHABET='abcdefghijklmnopqrstuvwxyz'

!! module parameters

     ! general codes
     !
     ! "\" character
     character(1), private, parameter :: C_ESCAPE=achar(27)
     !
     ! start ansi code, "\["
     character(2), private, parameter :: C_START=C_ESCAPE//'['
     !
     ! end ansi code, "m"
     character(1), private, parameter :: C_END='m'
     !
     ! clear all styles, "\[0m"
     character(4), private, parameter :: C_CLEAR=C_START//'0'//C_END

!!========================================================================
!!>>> declare local parameters                                         <<<
!!========================================================================

!! module parameters

     ! styles codes
     character(17), private, parameter :: STYLES(1:2,1:16)=reshape([&
         'BOLD_ON          ', '1  ', & ! Bold on.
         'ITALICS_ON       ', '3  ', & ! Italics on.
         'UNDERLINE_ON     ', '4  ', & ! Underline on.
         'INVERSE_ON       ', '7  ', & ! Inverse on: reverse foreground and background colors.
         'STRIKETHROUGH_ON ', '9  ', & ! Strikethrough on.
         'BOLD_OFF         ', '22 ', & ! Bold off.
         'ITALICS_OFF      ', '23 ', & ! Italics off.
         'UNDERLINE_OFF    ', '24 ', & ! Underline off.
         'INVERSE_OFF      ', '27 ', & ! Inverse off: reverse foreground and background colors.
         'STRIKETHROUGH_OFF', '29 ', & ! Strikethrough off.
         'FRAMED_ON        ', '51 ', & ! Framed on.
         'ENCIRCLED_ON     ', '52 ', & ! Encircled on.
         'OVERLINED_ON     ', '53 ', & ! Overlined on.
         'FRAMED_OFF       ', '54 ', & ! Framed off.
         'ENCIRCLED_OFF    ', '54 ', & ! Encircled off.
         'OVERLINED_OFF    ', '55 '  & ! Overlined off.
         ], [2,16])

!! module parameters

     ! colors codes: foreground colors
     character(15), private, parameter :: COLORS_FG(1:2,1:17)=reshape([&
         'BLACK          '  , '30 ', & ! black
         'RED            '  , '31 ', & ! red
         'GREEN          '  , '32 ', & ! green
         'YELLOW         '  , '33 ', & ! yellow
         'BLUE           '  , '34 ', & ! blue
         'MAGENTA        '  , '35 ', & ! magenta
         'CYAN           '  , '36 ', & ! cyan
         'WHITE          '  , '37 ', & ! white
         'DEFAULT        '  , '39 ', & ! default (white)
         'BLACK_INTENSE  '  , '90 ', & ! black intense
         'RED_INTENSE    '  , '91 ', & ! red intense
         'GREEN_INTENSE  '  , '92 ', & ! green intense
         'YELLOW_INTENSE '  , '93 ', & ! yellow intense
         'BLUE_INTENSE   '  , '94 ', & ! blue intense
         'MAGENTA_INTENSE'  , '95 ', & ! magenta intense
         'CYAN_INTENSE   '  , '96 ', & ! cyan intense
         'WHITE_INTENSE  '  , '97 '  & ! white intense
         ], [2,17])

     ! colors codes: background colors
     character(15), private, parameter :: COLORS_BG(1:2,1:17)=reshape([&
         'BLACK          '  , '40 ', & ! black
         'RED            '  , '41 ', & ! red
         'GREEN          '  , '42 ', & ! green
         'YELLOW         '  , '43 ', & ! yellow
         'BLUE           '  , '44 ', & ! blue
         'MAGENTA        '  , '45 ', & ! magenta
         'CYAN           '  , '46 ', & ! cyan
         'WHITE          '  , '47 ', & ! white
         'DEFAULT        '  , '49 ', & ! default (black)
         'BLACK_INTENSE  '  , '100', & ! black intense
         'RED_INTENSE    '  , '101', & ! red intense
         'GREEN_INTENSE  '  , '102', & ! green intense
         'YELLOW_INTENSE '  , '103', & ! yellow intense
         'BLUE_INTENSE   '  , '104', & ! blue intense
         'MAGENTA_INTENSE'  , '105', & ! magenta intense
         'CYAN_INTENSE   '  , '106', & ! cyan intense
         'WHITE_INTENSE  '  , '107'  & ! white intense
         ], [2,17])

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: pcs

     private :: color_index
     private :: style_index
     private :: upper

  contains ! encapsulated functionality

!!
!! @fun pcs
!!
!! colorize and stylize strings.
!!
  pure function pcs(string, fg, bg, style) result(cstr)
     implicit none

!! external arguments
     ! input string
     character(len=*), intent(in)           :: string

     ! foreground color definition
     character(len=*), intent(in), optional :: fg

     ! background color definition
     character(len=*), intent(in), optional :: bg

     ! style definition
     character(len=*), intent(in), optional :: style

!! local variables
     ! colorized string
     character(len=:), allocatable :: cstr

     ! counter
     integer(int32) :: i

!! [body

     cstr = string
     !
     ! if foreground color is specified
     if ( present(fg) ) then
         i = color_index(upper(fg))
         if ( i > 0 ) then
             cstr = C_START//trim(COLORS_FG(2,i))//C_END//cstr//C_CLEAR
         endif
     endif
     !
     ! if background color is specified
     if ( present(bg) ) then
         i = color_index(upper(bg))
         if ( i > 0 ) then
             cstr = C_START//trim(COLORS_BG(2,i))//C_END//cstr//C_CLEAR
         endif
     endif
     !
     ! if terminal style is specified
     if ( present(style) ) then
         i = style_index(upper(style))
         if ( i > 0 ) then
             cstr = C_START//trim(STYLES(2,i))//C_END//cstr//C_CLEAR
         endif
     endif

!! body]

     return
  end function pcs

!!
!! @fun color_index
!!
!! return the array-index corresponding to the queried color.
!!
!! because foreground and backround colors lists share the same name,
!! no matter what array is used to find the color index. thus, the
!! foreground array is used.
!!
  elemental function color_index(color)
     implicit none

!! external arguments
     ! color definition
     character(len=*), intent(in) :: color

!! local variables
     ! index into the colors arrays
     integer(int32) :: color_index

     ! counter
     integer(int32) :: c

!! [body

     color_index = 0
     do c=1, size(COLORS_FG, dim=2)
         if ( trim(COLORS_FG(1, c)) == trim(adjustl(color)) ) then
             color_index = c
             EXIT
         endif
     enddo

!! body]

     return
  end function color_index

!!
!! @fun style_index
!!
!! return the array-index corresponding to the queried style.
!!
  elemental function style_index(style)
     implicit none

!! external arguments
     ! style definition
     character(len=*), intent(in) :: style

!! local variables
     ! index into the styles array
     integer(int32) :: style_index

     ! counter
     integer(int32) :: s

!! [body

     style_index = 0
     do s=1, size(STYLES, dim=2)
         if ( trim(STYLES(1, s)) == trim(adjustl(style)) ) then
             style_index = s
             EXIT
         endif
     enddo

!! body]

     return
  end function style_index

!!
!! @fun upper
!!
!! return a string with all uppercase characters.
!!
  elemental function upper(string)
     implicit none

!! external arguments
     ! input string
     character(len=*), intent(in) :: string

!! local variables
     ! upper case string
     character(len=len(string)) :: upper

     ! characters counter
     integer :: n1
     integer :: n2

!! [body

     upper = string
     do n1=1, len(string)
         n2 = index(LOWER_ALPHABET, string(n1:n1))
         if ( n2 > 0 ) then
             upper(n1:n1) = UPPER_ALPHABET(n2:n2)
         endif
     enddo

!! body]

     return
  end function upper

  end module face
