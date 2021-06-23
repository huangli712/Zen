from h5 import *
from triqs.gf import *


h = HDFArchive("vasp.h5", "r")
gtau = h['dmft_output']['G_tau']['up_0']

for i in range(10001):
    print(i, " ", gtau.data[i,0,0].real, gtau.data[i,1,1].real, gtau.data[i,2,2].real, gtau.data[i,3,3].real, gtau.data[i,4,4].real)
