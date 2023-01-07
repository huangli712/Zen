/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "dmft.h"


DMFT::DMFT(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i):
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
		mdl.find_gfloc(seloc,gfloc);											if (mm) gfloc.write("gfloc", iter_cnt);
        update_phy();                                                           log("se_err_update");
		ImGreen g0effinv = gfloc.inverse() + seloc;						
		ImGreen g0eff = g0effinv.inverse();										if (mm) g0eff.write("g0eff", iter_cnt);
		ImGreen hb = find_hb(g0effinv);											if (mm) hb.write("hb", iter_cnt);
        ImGreen hb_c2 = rotate_hb(hb);											if (mm) hb_c2.write("hb_c2", iter_cnt);
		// if(!ifs || iter_cnt != 1) bth.number_bath_fit(hb_c2, iter_cnt);
		imp.update();
		norg_stage2.up_date_h0_to_solve(imp.h0);
		// norg_up_date_h0_with_solve_mutable(imp.h0, norg_stage1, norg_stage2);
		if(mm) std::cout << "The groundE_no-inter" << iofmt("sci") << norg_stage2.return_nointeractions_ground_state_energy(imp.h0) << std::endl;
        ImGreen g0imp(p.norbs * 2, p);	imp.find_g0(g0imp);						if (mm) g0imp.write("g0imp", iter_cnt);
		ImGreen gfimp(p.norbs * 2, p);norg_stage2.get_g_by_KCV_spup(gfimp);		if (mm) gfimp.write("gfimp", iter_cnt);// tested.
		// ImGreen gfimp_old(p.norbs * 2, p);norg_stage2.get_g_by_KCV(gfimp_old); 	if (mm) gfimp_old.write("gfimp_old", iter_cnt);// keep for test benchmark reason.
		// ImGreen gfimp_CF(2 * p.norbs, p);	norg_stage2.get_g_by_CF(gfimp_CF); 	if (mm) gfimp_CF.write("gfimp_CF", iter_cnt);// keep for test diagon benchmark reason.
        if(mm) std::cout << "The average of Green function" << iofmt("sci") << gfimp.particle_number() << std::endl;
        // norg_stage2.modify_by_krylov_correct_vec();
        ImGreen seimp(p.norbs * 2, p);
		seimp = g0imp.inverse() - gfimp.inverse();								if (mm) seimp.write("seimp", iter_cnt);
		append_se_err(seloc.error(seimp));
        seloc = seimp;
        // seloc = Cmplx(0.7) * seloc + Cmplx(0.3) * seimp;						if (mm) seloc.write("seloc", iter_cnt);
	}
	if (mm)	seloc.write("seloc");
	mdl.find_gfloc(seloc,gfloc);												if (mm) gfloc.write("gfloc");
    update_phy();                                                               log("se_err_converge");

    //ReGreen g0imp_re(p.norbs * 2, p);	imp.find_g0(g0imp_re);					if (mm) g0imp_re.write("g0imp_re");
    //ReGreen gimp_re(p.norbs * 2, p);	norg_stage2.find_g(gimp_re);			if (mm) gimp_re.write("gimp_re");
	//ReGreen se_re = g0imp_re.inverse() - gimp_re.inverse();					if (mm) se_re.write("se_re");
	//ReGreen gfloc_re(p.norbs * 2, p);	mdl.find_gfloc(se_re, gfloc_re);		if (mm) gfloc_re.write("gfloc_re");

}
// -------------------------------------------------- bethe lattice --------------------------------------------------


ImGreen DMFT::find_hb(const ImGreen& g0eff_inv) const
{
	ImGreen hb(g0eff_inv);
	MatCmplx h0loc_d = cmplx(mdl.hmlt0loc());
	for_Int(n, 0, p.num_omg) {
        hb[n] *= -1;
        hb[n] += dmat(hb.norbs, hb.z(n)) - h0loc_d;
    }
	//WRN(NAV(h0loc_d));
	return hb;
}

ImGreen DMFT::rotate_hb(const ImGreen &hb) const{
    ImGreen hb_r(hb);
    MatCmplx u = cmplx(direct_sum(p.c2u, p.c2u));
    for_Int(n, 0, hb.nomgs) {
        hb_r[n] = u.ct() * hb[n] * u;
    }
    return hb_r;
}


void DMFT::update_phy() {
    MatReal average = gfloc.particle_number();
    MatReal p_n_up = dmat(p.norbs,1.); 
    MatReal m_m_up = dmat(p.norbs, 1.);
    m_m_up[1][1] *= -1;
    m_m_up[2][2] *= -1;

    MatReal p_n_loc = direct_sum(p_n_up, -p_n_up);
    MatReal m_m_loc = direct_sum(m_m_up, m_m_up);

    p_n = 4.;
    m_m = 0.;
    for_Int(i, 0, p.norbs * 2) {
        p_n += p_n_loc[i][i] * average[i][i];
        m_m += m_m_loc[i][i] * average[i][i];
    }
    p_n /= 4;
    m_m /= 4;

    MatReal d_p_ud(p.norbs, 0.);	//ud=up-down
    MatReal p_p_ud(p.norbs, 0.);

    d_p_ud[0][1] = -1.;
    d_p_ud[1][0] = -1.;
    d_p_ud[0][2] = 1.;
    d_p_ud[2][0] = 1.;
    d_p_ud[1][3] = 1.;
    d_p_ud[3][1] = 1.;
    d_p_ud[2][3] = -1.;
    d_p_ud[3][2] = -1.;

    p_p_ud[0][1] = -1.;
    p_p_ud[1][0] = 1.;
    p_p_ud[0][2] = 1.;
    p_p_ud[2][0] = -1.;
    p_p_ud[1][3] = -1.;
    p_p_ud[3][1] = 1.;
    p_p_ud[2][3] = 1.;
    p_p_ud[3][2] = -1.;

	MatReal d_p_loc(2 * p.norbs, 2 * p.norbs, 0.);
	MatReal p_p_loc(2 * p.norbs, 2 * p.norbs, 0.);
    for_Int(i, 0, p.norbs) {
        for_Int(j, 0, p.norbs) {
            d_p_loc[i][j + p.norbs] = d_p_ud[i][j];
            d_p_loc[j + p.norbs][i] = d_p_ud.ct()[j][i];
            p_p_loc[i][j + p.norbs] = p_p_ud[i][j];
            p_p_loc[j + p.norbs][i] = p_p_ud.ct()[j][i];
        }
    }
    d_p = 0.;
    p_p = 0.;
    for_Int(row, 0, 2 * p.norbs) {
        for_Int(col, 0, 2 * p.norbs){
            d_p += average[row][col] * d_p_loc[col][row];
            p_p += average[row][col] * p_p_loc[col][row];
        }
    }
}