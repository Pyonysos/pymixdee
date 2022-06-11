"""
The pymixdee package is designed to help the scientist, statistician, etc., to construct appropriate mixture designs.
"""

import numpy as np
import pandas as pd
from typing import Union
from sympy.utilities.iterables import multiset_permutations



class MixD:
    """
    MixD class:
        attributes:
            nfact:              int, number of factors
            names:              str, list. names of the factors 
            row_limit:          int, maximum number of experiments
        methods:
            dirichlet:          
            simplex_centroid:
            scheffe_network:
    """

    def __init__(self, nfact : int=None, factor_names : list = None, row_limit : int = None):
        if (nfact is None) & (factor_names is None):
            raise ValueError('Missing arguments. Expecting value for either nfact or factor_names')
        self.nfact = nfact if nfact is not None else len(factor_names)
        self.names = factor_names
        self.row_limit = row_limit
        self.with_mixd = False
        self.with_constraints = False
        self.with_df_format = True


    def df_format(func):
        def inner(self, *args, **kwargs):
            mixd_array = func(self, *args, **kwargs)
            if self.with_df_format == False:
                return mixd_array
            elif self.names:
                return pd.DataFrame(mixd_array, columns=self.names)
            else:
                return mixd_array
        return inner

    @df_format    
    def dirichlet(self, size, alpha=None):
        """
        """
        if alpha == None:
            alpha = [1] * self.nfact
        elif 0 in alpha:
            alpha = [n if n > 0 else 1 for n in alpha]
        return np.random.default_rng().dirichlet(alpha, size)
    
    def add_constraints(self, constraints: Union[dict, list]):
        self.with_constraints = True

    def __add_center_points(self, mixd, ncenter):
        center_pt = np.ones((ncenter, self.nfact))/self.nfact
        return np.concatenate((mixd, center_pt), axis = 0)



    @df_format
    def simplex_centroid(self, ndegree=2, ncenter=1):
        """
        """

        ndegree = self.nfact - 1 if ndegree >= self.nfact else ndegree

        trim = np.tri(ndegree, self.nfact)
        trim /= np.sum(trim, axis=1).reshape(-1,1)

        for row in trim:
            permutations = np.array(list(multiset_permutations(row)))
            trim = np.concatenate((trim, permutations), axis = 0)

        if ncenter:
            mixd = self.__add_center_points(trim[ndegree:,:], ncenter)
        
        return mixd
        
    @df_format
    def scheffe_network(self, ndegree=2, ncenter=1):
        '''
        '''
        
        lattice = np.linspace(0, 1, ndegree, endpoint=False).reshape(-1,1)
        base = np.concatenate((lattice, 1-lattice, np.zeros((lattice.shape[0], self.nfact-2))), axis=1)

        for row in base:
            permutations = np.array(list(multiset_permutations(row)))
            base = np.concatenate((base, permutations), axis = 0)
        
        base = np.unique(base, axis=0)

        if ncenter > 0:
            mixd = self.__add_center_points(base, ncenter)
        
        return mixd
    
    
    '''
    in development
    '''
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
