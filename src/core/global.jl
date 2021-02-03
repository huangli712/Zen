#
# Project : pansy
# Source  : global.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : stable
# Comment :
#
# Last modified: 2021/02/03
#

"""
    I32 and I64

Alias of Integer type.
"""
const I32 = Int32
const I64 = Int64

"""
    F32 and F64

Alias of Float type.
"""
const F32 = Float32
const F64 = Float64

"""
    C32 and C64

Alias of Complex type.
"""
const C32 = ComplexF32
const C64 = ComplexF64

"""
    R32 and R64

Alias of Integer and Float types. Here `R` means Real.
"""
const R32 = Union{I32,F32}
const R64 = Union{I64,F64}

"""
    N32 and N64

Alias of Integer, Float, and Complex types. Here `N` means Number.
"""
const N32 = Union{I32,F32,C32}
const N64 = Union{I64,F64,C64}

"""
    __LIBNAME__

Name of this julia package.
"""
const __LIBNAME__ = "ZEN Framework"

"""
    __VERSION__

Version of this julia package.
"""
const __VERSION__ = "0.0.8@devel"

"""
    __RELEASE__

Release date of this julia package.
"""
const __RELEASE__ = "2021/02"

"""
    __AUTHORS__

Core authors of this julia package.
"""

#
# Remarks:
#
# The Array's element should be a NamedTuple object, such as:
#     (name = "name", email = "email").
#
const __AUTHORS__ = [(name = "Li Huang", email = "lihuang.dmft@gmail.com")]
