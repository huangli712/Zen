/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2023
*/

#include "apizen.h"


APIzen::APIzen(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i):
	mm(mm_i), p(prmtr_i), mdl(mdl_i), iter_cnt(0), gfloc(p.norbs * 2, p),
	se_err({ 1.e99,1.e98,1.e97 ,1.e96 ,1.e95}), seloc(p.norbs * 2, p),
	//se_err({ 1.e99 }),
	p_n(0.),
	m_m(0.),
	d_p(0.),
	p_p(0.)
{
	// make random seed output together
	{ mm.barrier(); SLEEP(1); }
	Bath bth(mm, p);
	// IFS ifs("bi/001.fvb.txt");	if (ifs) bth.fvb = bth.bth_read_fvb("bi/001.fvb.txt");
	Impurity imp(mm, p, bth, mdl);
	// NORG norg_stage1(mm, p); shift_norg_restrain(); 
	NORG norg_stage2(mm, p);
	// seloc = ImGreen(p.norbs * 2, p, "bi/input_seloc.txt");									if (mm) seloc.write("seloc", iter_cnt);

	while (iter_cnt < p.iter_max && !converged()) 
	{
		++iter_cnt;
		// mdl.find_gfloc(seloc,gfloc);												if (mm) gfloc.write("gfloc", iter_cnt);
        // update_phy();                                                           	log("se_err_update");
		// ImGreen g0effinv = gfloc.inverse() + seloc;						
		// ImGreen g0eff = g0effinv.inverse();										if (mm) g0eff.write("g0eff", iter_cnt);
		// ImGreen hb = find_hb(g0effinv);											if (mm) hb.write("hb", iter_cnt);
        // ImGreen hb_c2 = rotate_hb(hb);											if (mm) hb_c2.write("hb_c2", iter_cnt);
		// // if(!ifs || iter_cnt != 1) bth.number_bath_fit(hb_c2, iter_cnt);
		// imp.update();
		// norg_stage2.up_date_h0_to_solve(imp.h0);
		// // norg_up_date_h0_with_solve_mutable(imp.h0, norg_stage1, norg_stage2);
		// if(mm) std::cout << "The groundE_no-inter" << iofmt("sci") << norg_stage2.return_nointeractions_ground_state_energy(imp.h0) << std::endl;
        // ImGreen g0imp(p.norbs * 2, p);	imp.find_g0(g0imp);						if (mm) g0imp.write("g0imp", iter_cnt);
		// ImGreen gfimp(p.norbs * 2, p);norg_stage2.get_g_by_KCV_spup(gfimp);		if (mm) gfimp.write("gfimp", iter_cnt);// tested.
		// // ImGreen gfimp_old(p.norbs * 2, p);norg_stage2.get_g_by_KCV(gfimp_old); 	if (mm) gfimp_old.write("gfimp_old", iter_cnt);// keep for test benchmark reason.
		// // ImGreen gfimp_CF(2 * p.norbs, p);	norg_stage2.get_g_by_CF(gfimp_CF); 	if (mm) gfimp_CF.write("gfimp_CF", iter_cnt);// keep for test diagon benchmark reason.
        // if(mm) std::cout << "The average of Green function" << iofmt("sci") << gfimp.particle_number() << std::endl;
        // // norg_stage2.modify_by_krylov_correct_vec();
        // ImGreen seimp(p.norbs * 2, p);
		// seimp = g0imp.inverse() - gfimp.inverse();								if (mm) seimp.write("seimp", iter_cnt);
		// append_se_err(seloc.error(seimp));
        // seloc = seimp;
        // seloc = Cmplx(0.7) * seloc + Cmplx(0.3) * seimp;						if (mm) seloc.write("seloc", iter_cnt);
	}
	if (mm)	seloc.write("seloc");
	mdl.find_gfloc(seloc,gfloc);												if (mm) gfloc.write("gfloc");
    // update_phy();                                                               log("se_err_converge");

    //ReGreen g0imp_re(p.norbs * 2, p);	imp.find_g0(g0imp_re);					if (mm) g0imp_re.write("g0imp_re");
    //ReGreen gimp_re(p.norbs * 2, p);	norg_stage2.find_g(gimp_re);			if (mm) gimp_re.write("gimp_re");
	//ReGreen se_re = g0imp_re.inverse() - gimp_re.inverse();					if (mm) se_re.write("se_re");
	//ReGreen gfloc_re(p.norbs * 2, p);	mdl.find_gfloc(se_re, gfloc_re);		if (mm) gfloc_re.write("gfloc_re");

}
// ----------------------------------------------------------------------------------------------------
