It deals with the error and exception.

## Type

subroutines

## Source

`src/s_error.f90`

## Usage

**(1)** Display the error, exception, or message on the terminal.

```fortran
subroutine s_print_error(sub, msg)
subroutine s_print_exception(sub, msg)
subroutine s_print_message(sub, msg)
```

where `sub` is the name of subroutine or function that emits the error, exception, or message, `msg` explains the error, exception, or message.

**(2)** Write the error, exception, or message to external file.

```fortran
subroutine s_write_error(sub, msg, file_unit)
subroutine s_write_exception(sub, msg, file_unit)
subroutine s_write_message(sub, msg, file_unit)
```

where `sub` is the name of subroutine or function that emits the error, exception, or message, `msg` explains the error, exception, or message, `file_unit` denotes the handler of external file.
