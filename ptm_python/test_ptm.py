import unittest
import numpy as np
import ptm_postprocessing as ppo


__all__ = ['postProcTests']


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


if __name__ == "__main__":
    unittest.main()
