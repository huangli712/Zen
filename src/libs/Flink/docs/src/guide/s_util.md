It provides some utility subroutines and functions.

## Type

subroutines

## Source

`src/s_util.f90`

## Usage

**(1)** C-like assertion.

```fortran
subroutine s_assert(condition)
subroutine s_assert2(condition, message)
```

where `condition` is the logical condition that we have to assert, `message` is the additional message.

**(2)** Sort algorithm.

```fortran
subroutine s_sorter1_i(nsize, list)
subroutine s_sorter1_d(nsize, list)
subroutine s_sorter2_i(nsize, list, indx)
subroutine s_sorter2_d(nsize, list, indx)
subroutine s_sorter3_i(nsize, list)
subroutine s_sorter3_d(nsize, list)
subroutine s_sorter4_i(nsize, list, indx)
subroutine s_sorter4_d(nsize, list, indx)
subroutine s_qsorter(nsize, list)
subroutine s_qscorer(pstart, pend, nsize, list)
```

The `s_sorter1()` and `s_sorter2()` subroutines implement the bubble algorithm. The `s_sorter3()` and `s_sorter4()` subroutines implement the heap algorithm. And the `s_qsorter()` subroutine implements the quick sort algorithm. The `s_qscorer()` subroutine is called by the `s_qsorter()` internally. DO NOT call it directly!

Here, `nsize` is size of the list, `list` is the dataset to be sorted. In `s_sorter2()` and `s_sorter4()`, the index of the dataset (`indx`) is also sorted.

!!! note

    `_i` means integer version, `_d` real(dp) version.

**(3)** Combination algrbra.

```fortran
subroutine s_combination(ntiny, nlarg, value)
```

```fortran
     ! the small number
     integer, intent(in)  :: ntiny

     ! the large number
     integer, intent(in)  :: nlarg

     ! result value of the combination algebra
     integer, intent(out) :: value
```

It calculates combination algebra.

**(4)** String manipulation.

```fortran
subroutine s_str_upcase(s)
subroutine s_str_lowcase(s)
```

```fortran
     ! input/output string
     character(len=*), intent(inout) :: s
```

They are used to convert the input string `s` to upcase or lowcase.

```fortran
subroutine s_str_count(string, substr, count)
```

```fortran
     ! string to be examined
     character(len=*), intent(in) :: string

     ! substring in question
     character(len=*), intent(in) :: substr

     ! return value, number of occurrences
     integer, intent(out) :: count
```

It return the number of times a substring occurs.

```fortran
subroutine s_str_compress(string)
```

```fortran
     ! character string to be compressed.
     character(len=*), intent(inout) :: string
```

It return a copy of an input string with all whitespace (spaces and tabs) is removed.

**(5)** Date time manipulation.

```fortran
subroutine s_time_builder(date_time_string)
```

```fortran
     ! output date and time
     character(len=20), intent(out) :: date_time_string
```

It returns a string containing date and time in human readable format.

```fortran
subroutine s_time_analyzer(time_iter, time_niter)
```

```fortran
     ! time used in this iteration
     real(dp), intent(in) :: time_iter

     ! time used in total iteration
     real(dp), intent(in) :: time_niter
```

It is used to print the iteration timing information about continuous time quantum Monte Carlo quantum impurity solver.
