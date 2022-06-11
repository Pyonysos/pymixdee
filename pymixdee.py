"""
The pymixdee package is designed to help the scientist, statistician, etc., to construct appropriate mixture designs.
"""

import numpy as np
from typing import Union
from sympy.utilities.iterables import multiset_permutations


class MixD:
    """
    MixD class:
        attributes:
            nfact:      int, number of factors
            names:      str, list. names of the factors 
            row_limit:  int, maximum number of experiments
        methods:
    """

    def __init__(self, nfact : int, factor_names : list = None, row_limit : int = None):
        self.nfact = nfact
        self.names = factor_names
        self.row_limit = row_limit
        self.with_mixd = False
        self.with_constraints = False

    def __diag(self, n):
        """
        """
        return np.eye(n)
    
    def __dirichlet(self, size, alpha=None):
        """
        """
        if alpha == None:
            alpha = [1] * self.nfact
        return np.random.default_rng().dirichlet(alpha, size)
    
    def __combinations(self, nfact):
        trim = np.tri(nfact, nfact, 1)
        trim /= trim.sum(axis=1)
        for row in trim:
            trim = np.vstack(trim, np.array(multiset_permutations(row)))
        return trim
    
    def add_constraints(self, constraints: Union[dict, list]):
        self.with_constraints = True
    
    def simplex_centroid(self, nfact, domain=None):
        """
        """
        

    def scheffe_network(self):
        """
        """
        ...
    
    def generate(self, design = None):
        """
        """
        if design in ('centroid', 'sc'):
            simplex_centroid()

        elif design in ('lattice', 'scheffe', 'sn'):
            scheffe_network()

        elif design in ('doptimal'):
            pass

        elif design in ('dirichlet'):
            pass

    def d_optimal(self):
        """
        """
        ...


    def plot_experiments(self):
        """
        """
        if self.with_mixd:
            fig = ff.create_ternary_contour(x, y, pole_labels=self.names, interp_mode='cartesian', colorscale='Viridis', showscale=True, ncontours=ncontours)
            fig.show()
