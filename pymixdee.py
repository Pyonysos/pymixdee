import numpy as np

class MixD:
    def __init__(self, nfact=None, factor_names=None):
        self.nfact = nfact
        self.names = factor_names
    
    def __diag(self, n):
        return np.eye(n)
    
    def __dirichlet(self, size, alpha=None):
        if alpha == None:
            alpha = [1] * self.nfact
        return np.random.default_rng().dirichlet(alpha, size)
    
    def simplex_centroid(self, nfact, domain=None):
        ...

    def scheffe_network(self):
        ...
    
    def mixD(self, type=None):
        if type == 'centroid':
            simplex_centroid()

        if type == 'lattice':
            scheffe_network()

    def d_optimal(self):
        ...

    def cmixd(self, nfact, ):
        ...



    def plot(self):
        ...
