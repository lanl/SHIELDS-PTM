import unittest
import numpy as np
import ptm_postprocessing as ppo
import ptm_tools as ptt


__all__ = ['postProcTests', 'toolsTests']


class postProcTests(unittest.TestCase):
    def setUp(self):
        super(postProcTests)
        self.pp = ppo.ptm_postprocessor()

    def tearDown(self):
        super(postProcTests)

    def test_dist_max(self):
        res = self.pp.test_distributions(verbose=False)['Maxwell']
        np.testing.assert_almost_equal(res[0], res[1])

    def test_dist_int_max(self):
        """Integral of Maxwell should equal 1"""
        res = self.pp.test_distributions(verbose=False)['Maxwell']
        np.testing.assert_almost_equal(1, res[2])

    def test_dist_jut(self):
        res = self.pp.test_distributions(verbose=False)['Juttner']
        np.testing.assert_almost_equal(res[0], res[1])

    def test_dist_int_jut(self):
        """Integral of Maxwell-Juttner should equal 1"""
        res = self.pp.test_distributions(verbose=False)['Juttner']
        np.testing.assert_almost_equal(1, res[2])

    def test_dist_kap(self):
        res = self.pp.test_distributions(verbose=False)['Kappa']
        np.testing.assert_almost_equal(res[0], res[1])

    def test_dist_int_kap(self):
        """Integral of Kappa should equal 1"""
        res = self.pp.test_distributions(verbose=False)['Kappa']
        np.testing.assert_almost_equal(1, res[2])

    def test_dist_rkt(self):
        res = self.pp.test_distributions(verbose=False)['relKT']
        np.testing.assert_almost_equal(res[0], res[1])

    def test_dist_int_rkt(self):
        """Integral of relativistic Kappa-type should equal 1"""
        res = self.pp.test_distributions(verbose=False)['relKT']
        np.testing.assert_almost_equal(1, res[2])


class toolsTests(unittest.TestCase):
    def setUp(self):
        super()
        self.vstorm_2010 = ptt.StormerVertical()
        self.cstorm_2010 = ptt.StormerCutoff()

    def tearDown(self):
        super()

    def test_proton_roundtrip(self):
        '''cutoff rigidity for given energy should equal energy for given rigidity'''
        pred_e = []
        init_e = [1, 10, 100, 1000]
        for ie in init_e:
            rig = ptt.Proton(ie).getRigidity()  # 100 MeV
            prot = ptt.Proton.fromRigidity(rig)
            pred_e.append(prot.energy)
        np.testing.assert_array_almost_equal(pred_e, init_e)

    def test_Stormer_gen_to_vert(self):
        '''Stormer cutoff with zenith angle zero should be same as Stormer vertical approx'''
        for l_val in [2.5, 4, 6.6]:
            v_rig = self.vstorm_2010.cutoff_at_L(l_val)
            c_rig = self.cstorm_2010.cutoff_at_L(l_val, 0, 0)
            self.assertAlmostEqual(v_rig, c_rig)

    def test_Stormer_azi_at_z0(self):
        '''Stormer cutoff should not depend on azimuth if zenith angle is zero'''
        l_val = 5
        init_c = [self.cstorm_2010.cutoff_at_L(l_val, 0, 0)] * 4
        esti_c = [self.cstorm_2010.cutoff_at_L(l_val, 0, az) for az in (0, 45, 180, 199.3)]
        np.testing.assert_array_almost_equal(esti_c, init_c)

    def test_Stormer_zenith_nadir(self):
        '''Stormer cutoff should not be higher for nadir direction than zenith'''
        l_val = 5
        zen_c = self.cstorm_2010.cutoff_at_L(l_val, 0, 90)
        nad_c = self.cstorm_2010.cutoff_at_L(l_val, 180, 90)
        self.assertLessEqual(nad_c, zen_c)

    def test_Stormer_east_west(self):
        '''Stormer cutoff should be lower for west than east'''
        l_val = 5
        c_east = self.cstorm_2010.cutoff_at_L(l_val, 90, 90)
        c_west = self.cstorm_2010.cutoff_at_L(l_val, 90, 270)
        self.assertLessEqual(c_west, c_east)

    def test_Stormer_east_west_array(self):
        '''Stormer cutoff should be lower for west than east and work with array input'''
        l_val = [5, 6]
        c_east = self.cstorm_2010.cutoff_at_L(l_val, 90, 90, as_energy=True)
        c_west = self.cstorm_2010.cutoff_at_L(l_val, 90, 270, as_energy=True)
        np.testing.assert_array_less(c_west, c_east)

    def test_Stormer_all_array(self):
        '''Stormer cutoff work with array/list input for all 3 vars'''
        l_val = [5, 6]
        zen = [90, 90]
        azi = [90, 270]
        c_east = self.cstorm_2010.cutoff_at_L(l_val[0], 90, 90, as_energy=True)
        c_west = self.cstorm_2010.cutoff_at_L(l_val[1], 90, 270, as_energy=True)
        indiv = np.array([c_east, c_west])
        togeth = self.cstorm_2010.cutoff_at_L(l_val, zen, azi, as_energy=True)
        np.testing.assert_array_almost_equal(indiv, togeth)

    def test_Stormer_all_array(self):
        '''Stormer cutoff work with array/list input for all 3 vars'''
        l_val = [5, 6, 7]
        zen = 90
        azi = 270
        c1 = self.cstorm_2010.cutoff_at_L(l_val[0], zen, azi, as_energy=True)
        c2 = self.cstorm_2010.cutoff_at_L(l_val[1], zen, azi, as_energy=True)
        c3 = self.cstorm_2010.cutoff_at_L(l_val[2], zen, azi, as_energy=True)
        indiv = np.array([c1, c2, c3])
        togeth = self.cstorm_2010.cutoff_at_L(l_val, zen, azi, as_energy=True)
        np.testing.assert_array_almost_equal(indiv, togeth)

    def test_l_to_lati_roundtrip(self):
        '''dipole L to invariant lat should be reversible'''
        l_val = [1, 2, 4, 8]
        lati = [ptt.invariant_latitude_from_l(ll) for ll in l_val]
        pred_l = [ptt.l_from_invariant_latitude(il) for il in lati]
        np.testing.assert_array_almost_equal(l_val, pred_l)

    def test_l_to_lati_roundtrip_array(self):
        '''dipole L to invariant lat should be reversible'''
        l_val = [1, 2, 4, 8]
        lati = ptt.invariant_latitude_from_l(l_val)
        pred_l = ptt.l_from_invariant_latitude(lati)
        np.testing.assert_array_almost_equal(l_val, pred_l)

    def test_Stormer_vertical_cutoff_list(self):
        lvals = [2.5, 4, 6.6]
        v_rig_1 = [self.vstorm_2010.cutoff_at_L(l_val, as_energy=True) for l_val in lvals]
        v_rig_2 = self.vstorm_2010.cutoff_at_L(lvals, as_energy=True)
        np.testing.assert_array_almost_equal(v_rig_1, v_rig_2)

if __name__ == "__main__":
    unittest.main()
