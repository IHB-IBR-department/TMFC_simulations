import numpy as np
from neurolib.models.wc import WCModel


def simple_one_node(exc_ext=1.5):
    wc = WCModel()
    wc.params['duration'] = 2.0 * 1000
    wc.params['exc_ext'] = exc_ext
    wc.params['c_inhexc'] = 12
    wc.run()
    return wc


def ito_version(K_gl=4, exc_ext=-5):
    A = np.array([[0, 1], [1, 0]])
    D = np.array([[0, 0], [0, 0]])
    wc = WCModel(Cmat=A, Dmat=D)
    wc.params['duration'] = 1000
    wc.params['exc_ext'] = exc_ext
    wc.params['K_gl'] = K_gl
    wc.params['dt'] = 10
    wc.params['c_inhexc'] = 0
    wc.params['sigma_ou'] = 0
    wc.params['mu_exc'] = 2
    wc.params['a_exc'] = 1
    wc.params['c_excexc'] = 2
    wc.params['tau_exc'] = 1
    wc.run()
    return wc
