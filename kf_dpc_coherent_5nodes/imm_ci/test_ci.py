import unittest
import numpy as np
from fusion import forecast, fuse, kernel
from imm_filter import fit, predict


class Checks(unittest.TestCase):
    def setUp(self):
        self.W=np.array([[.5,.5,0],[.5,0,.5],[0,.5,.5]])
        self.C=kernel(self.W,4)

    def test_identical_input_is_not_counted_multiple_times(self):
        x=np.broadcast_to(np.arange(6.),(2,3,6)).copy()
        P=np.broadcast_to(np.diag(np.arange(1,7.)),(2,3,6,6)).copy()
        y,Q=fuse(x,P,self.C,'ci')
        np.testing.assert_allclose(y,x,atol=1e-13)
        np.testing.assert_allclose(Q,P,atol=1e-13)

    def test_equal_covariance_reduces_to_state_average(self):
        x=np.arange(36.).reshape(2,3,6)
        P=np.broadcast_to(np.eye(6),(2,3,6,6))
        y,Q=fuse(x,P,self.C,'ci')
        np.testing.assert_allclose(y,np.einsum('ij,tjd->tid',self.C,x),atol=1e-13)
        np.testing.assert_allclose(Q,P,atol=1e-13)

    def test_finite_round_matches_explicit_neighbour_information_exchange(self):
        rng=np.random.default_rng(1);A=rng.normal(size=(2,3,6,6))
        P=A@A.swapaxes(-1,-2)+np.eye(6)
        x=rng.normal(size=(2,3,6));Y=np.linalg.inv(P);y=np.einsum('tnij,tnj->tni',Y,x)
        for _ in range(4):
            Y=np.einsum('ij,tjab->tiab',self.W,Y);y=np.einsum('ij,tjd->tid',self.W,y)
        xf,Pf=fuse(x,P,self.C,'ci')
        np.testing.assert_allclose(Pf,np.linalg.inv(Y),atol=1e-13)
        np.testing.assert_allclose(xf,np.linalg.solve(Y,y[...,None])[...,0],atol=1e-13)
        self.assertGreater(np.linalg.eigvalsh(Pf).min(),0)

    def test_forecast_zero_and_nonzero_mean_matches_existing_imm(self):
        t=np.arange(15)*.2;z=np.zeros((15,3,3));z[...,0]=t[:,None]*2
        R=np.broadcast_to(np.eye(3),(15,3,3,3))
        h=fit(t,z,R,retain_branch_covariances=True)
        x,P=forecast(h,0)
        np.testing.assert_allclose(x,h['state'][...,:6],atol=1e-13)
        np.testing.assert_allclose(P,h['covariance'][...,:6,:6],atol=1e-13)
        x,P=forecast(h,.177)
        np.testing.assert_allclose(x,predict(h,.177),atol=1e-13)
        self.assertGreater(np.linalg.eigvalsh(P).min(),0)

if __name__=='__main__':unittest.main()
