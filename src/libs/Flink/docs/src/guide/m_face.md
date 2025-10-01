It supports terminal escape sequences to produce styled and colored terminal output.

## Type

module

## Source

`src/m_face.f90`

## Usage

**(1)** Colors

The foreground and background colors could be:

* black
* red
* green
* yellow
* blue
* magenta
* cyan
* white
* default (it is actually black)

**(2)** Styles

The possible styles could be:

* bold_on
* italics_on
* underline_on
* inverse_on
* strikethrough_on
* bold_off
* italics_off
* underline_off
* inverse_off
* strikethrough_off
* framed_on
* encircled_on
* overlined_on
* framed_off
* encircled_off
* overlined_off

**(3)** Function

```fortran
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
```

The possible values for the `fg`, `bg`, and `style` arguments have been presented above.

**(4)** Sample codes

```fortran
program test
    use face

    implicit none

    print *, pcs("hello world", 'white', 'black', 'bold_on')
end program test
```
