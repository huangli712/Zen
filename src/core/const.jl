#
# project : pansy
# source  : const.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : stable
# comment :
#
# last modified: 2020/12/29
#

"""
    I32 and I64

Alias of integer type
"""
const I32 = Int32
const I64 = Int64

"""
    F32 and F64

Alias of float type
"""
const F32 = Float32
const F64 = Float64

"""
    C32 and C64

Alias of complex type
"""
const C32 = ComplexF32
const C64 = ComplexF64

"""
    __libname__

Name of this julia package
"""
const __libname__ = "ZEN"

"""
    __version__

Version of this julia package
"""
const __version__ = "0.0.5@d"

"""
    __release__

Release date of this julia package
"""
const __release__ = "2020/12"

"""
    __authors__

Core authors of this julia package
"""
const __authors__ = [("Li Huang", "lihuang.dmft@gmail.com")]
