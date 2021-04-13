from h5 import *

f = open("vasp.h5.txt", "a")

h = HDFArchive("vasp.h5", "r")

for i in h['dft_input'].items():
    print(i, file = f)

f.close()
