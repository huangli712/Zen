import triqs_dft_tools.converters.plovasp.converter as plo_converter
from triqs_dft_tools.converters.vasp import *
from triqs_dft_tools.sumk_dft_tools import SumkDFTTools

plo_converter.generate_and_output_as_text('plo.cfg', vasp_dir='./')

# create Converter
Converter = VaspConverter(filename = 'vasp')

# run the converter
Converter.convert_dft_input()

SK = SumkDFTTools(hdf_file='vasp.h5', use_dft_blocks = False)

Sigma = SK.block_structure.create_gf(beta=40)
SK.put_Sigma([Sigma])
G = SK.extract_G_loc()
SK.analyse_block_structure_from_gf(G, threshold = 1e-3)
for i_sh in range(len(SK.deg_shells)):
    num_block_deg_orbs = len(SK.deg_shells[i_sh])
    mpi.report('found {0:d} blocks of degenerate orbitals in shell {1:d}'.format(num_block_deg_orbs, i_sh))
    for iblock in range(num_block_deg_orbs):
        mpi.report('block {0:d} consists of orbitals:'.format(iblock))
        for keys in list(SK.deg_shells[i_sh][iblock].keys()):
            mpi.report('  '+keys)
print("hehe")
