import unittest
import numpy as np
from imm_filter import transition, jacobian, fit, predict, Parameters, mix


class Checks(unittest.TestCase):
    def test_ct_known_circle_and_zero_limit(self):
        x=np.array([[0.,0.,10.,15.,0.,0.,.075]])
        y=transition(x,2.,2)
        np.testing.assert_allclose(y[0,:2],[200*np.sin(.15),200*(1-np.cos(.15))])
        x[0,6]=0
        np.testing.assert_allclose(transition(x,.2,2)[...,:6],transition(x,.2,1)[...,:6])
        self.assertAlmostEqual(jacobian(x,.2,2)[0,1,6],.5*15*.2**2)

    def test_mixture_between_model_uncertainty(self):
        x=np.zeros((1,2,7));x[0,1,0]=10
        P=np.broadcast_to(np.eye(7),(1,2,7,7)).copy()
        m,c,_=mix(x,P,np.array([[.5,.5]]),np.full((2,2),.5))
        np.testing.assert_allclose(m[...,0],5)
        np.testing.assert_allclose(c[...,0,0],26)

    def test_causal_covariance_and_probabilities(self):
        t=np.arange(51)*.2
        z=np.broadcast_to((np.array([3.,-2.,1.])+t[:,None]*[4.,2.,.2])[:,None],(51,2,3)).copy()
        R=np.broadcast_to(np.eye(3)*.01,(51,2,3,3))
        h=fit(t,z,R); prefix=fit(t[:20],z[:20],R[:20])
        np.testing.assert_array_equal(h['state'][:20],prefix['state'])
        np.testing.assert_allclose(h['probabilities'].sum(-1),1,atol=1e-13)
        self.assertGreater(np.linalg.eigvalsh(h['covariance']).min(),0)
        np.testing.assert_allclose(h['state'][-1,0,3:6],[4,2,.2],atol=.05)
        np.testing.assert_allclose(predict(h,0),h['state'][...,:6],atol=1e-13)

    def test_stationary_is_not_forced_to_origin(self):
        t=np.arange(151)*.2; rng=np.random.default_rng(8)
        z=np.array([120.,-80.,50.])+rng.normal(0,.1,(151,1,3))
        h=fit(t,z,np.broadcast_to(np.eye(3)*.01,(151,1,3,3)))
        np.testing.assert_allclose(h['state'][-1,0,:3],[120,-80,50],atol=.1)
        self.assertGreater(h['probabilities'][-30:,:,0].mean(),.5)

    def test_unknown_turn_rate_learned_from_positions(self):
        t=np.arange(151)*.2; w=.075
        z=np.stack([200*np.sin(w*t),-200*np.cos(w*t),np.full(len(t),120.)],-1)[:,None]
        h=fit(t,z,np.broadcast_to(np.eye(3)*.0001,(len(t),1,3,3)),enabled=(2,))
        self.assertEqual(h['state'][0,0,6],0.)
        self.assertAlmostEqual(h['state'][-1,0,6],w,places=5)
        np.testing.assert_allclose(h['state'][-1,0,3:6],[15*np.cos(w*t[-1]),15*np.sin(w*t[-1]),0],atol=.001)

if __name__=='__main__': unittest.main()
