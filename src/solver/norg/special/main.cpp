
#include "specs.h"
#include "prmtr.h"
#include "model.h"
#include "apizen.h"

int main(int argc, char* argv[])
{
	using namespace std;
	MPI_Init(&argc, &argv);
	MyMpi mm;
	if (mm) std::cout << "\n\n Version: v0.2.13 @ 2023.02.22 (running "<< present() <<")\n\n" << std::endl;
	if (mm) cout << NAV(pwd()) << endl; 
	use_mkl(mm);
	// if (mm) io_init();
	clock_t t_program_bgn;
	if (mm) TIME_BGN("program", t_program_bgn);
	Prmtr prmtr(mm);

	// Model mdl(mm,prmtr);
	{// For test NORG only.
		// Bath bth(mm, prmtr);
		// Impurity imp(mm, prmtr, bth, mdl);
		// NORG norg(mm,prmtr);
		// ImGreen g0(2 * prmtr.norbs, prmtr);	imp.find_g0(g0); if (mm) g0.write("Theor_g0", 1);
		// norg.up_date_h0_to_solve(imp.h0);
		
		// ImGreen norg_KCV(2 * prmtr.norbs, prmtr);	norg.get_g_by_KCV(norg_KCV);	if (mm) norg_KCV.write("NORG_KCV", 2);
		// // ImGreen norg_CF(2 * prmtr.norbs, prmtr); 	norg.get_g_by_CF(norg_CF); 		if (mm) norg_CF.write("NORG_CF", 3);

		// // if(mm) WRN(NAV(norg.return_nointeractions_ground_state_energy(imp.h0)));
		// if(mm) std::cout << "The groundE_pre" << iofmt("sci") << norg.return_nointeractions_ground_state_energy(imp.h0) << std::endl;
	}
	// DMFT dmft(mm, prmtr, mdl);
    // Phy phy(mm, prmtr, mdl);
	APIzen norg(mm, prmtr, "solver", 0);
	
	// ImGreen gfimp(prmtr.nband, prmtr, "gfimp.txt");		if(mm) WRN(NAV(gfimp.particle_number().diagonal()));
	// ImGreen g0imp(prmtr.nband, prmtr, "gfimp0.txt");	if(mm) WRN(NAV(g0imp.particle_number().diagonal()));
	// ImGreen seimp(prmtr.nband, prmtr);	seimp = g0imp.inverse() - gfimp.inverse();	if (mm) seimp.write_zen("seimp");
    if (mm)	TIME_END("program", t_program_bgn);
    MPI_Finalize();
	return 0;
}