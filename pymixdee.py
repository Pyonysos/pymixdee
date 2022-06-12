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
        self.constraints = {}
        
        self.with_df_format = True

    @df_format
    def shuffle(self, mixd):
        mixd = np.array(mixd)
        np.random.default_rng().shuffle(mixd)
        return mixd

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
    def dirichlet(self, size: int, alpha: list[int]=None):
        """
        generate a series of mixture following a dirichlet distribution
        size:   int, number a mix to generate
        alpha:  list of int
        """
        if alpha == None:
            alpha = [1] * self.nfact
        elif 0 in alpha:
            alpha = [n if n > 0 else 1 for n in alpha]
        return np.random.default_rng().dirichlet(alpha, size)
    
    def add_lower_constraints(self, mixd, constraints : list):
        '''
        add lower constraints to the factors the mix design
          mixd: pd.Dataframe or np.array. mix plan
          constraints: list of float. 
        '''
        constraints = np.array(constraints).reshape(1, -1)
        mixd *= (1 - constraints.sum())
        mixd += constraints
        return mixd

    def __add_center_points(self, mixd, ncenter):
        '''
        add center points to mixture design
        mixd:     np.array or pd.DataFrame, mixture design
        ncenter:  int, number of central points to be added
        '''
        center_pt = np.ones((ncenter, self.nfact))/np.array(self.nfact)
        return np.concatenate((mixd, center_pt), axis = 0)

    @df_format
    def simplex_centroid(self, ndegree=2, ncenter=1, lower: list= None):
        """
        generate centroid simplex
        ndegree:    int, default 2, number of degree of design
        ncenter:    int, default 1, number of central points to be added
        lower:      list of float, default None, lower constraints to each factor
        """
        ndegree = self.nfact - 1 if ndegree >= self.nfact else ndegree

        trim = np.tri(ndegree, self.nfact)
        trim /= np.sum(trim, axis=1).reshape(-1,1)

        for row in trim:
            permutations = np.array(list(multiset_permutations(row)))
            trim = np.concatenate((trim, permutations), axis = 0)

        if ncenter:
            mixd = self.__add_center_points(trim[ndegree:,:], ncenter)
        print('number of experiments =' , mixd.shape[0])
    
        if lower is not None:
            mixd = self.add_lower_constraints(mixd, lower)

        return mixd
        
    @df_format
    def scheffe_network(self, ndegree=2, ncenter=1, lower: list=None):
        """
        generate scheffe network
        ndegree:    int, default 2, number of degree of design
        ncenter:    int, default 1, number of central points to be added
        lower:      list of float, default None, lower constraints to each
        """
        
        lattice = np.linspace(0, 1, ndegree, endpoint=False).reshape(-1,1)
        base = np.concatenate((lattice, 1-lattice, np.zeros((lattice.shape[0], self.nfact-2))), axis=1)

        for row in base:
            permutations = np.array(list(multiset_permutations(row)))
            base = np.concatenate((base, permutations), axis = 0)

        duplicates = []

        for n in range(base.shape[0]):
            for m in range(n+1, base.shape[0]):
                if np.all(np.isclose(base[n], base[m])):
                    duplicates.append(m)
        base = base[[n for n in range(base.shape[0]) if n not in duplicates]]

        #base = base[indexes]

        if ncenter > 0:
            mixd = self.__add_center_points(base, ncenter, dtype=dtype)
        print('number of experiments =' , mixd.shape[0])
        
        
        if lower is not None:
            mixd = self.add_lower_constraints(mixd, lower)
        return mixd
    
    def export(self, mixd:pd.DataFrame, filename:str, extension: str='xlsx')
        if extension in ('xlsx', 'excel'):
          self.export_to_excel(mixd, filename)
        elif extension == 'csv':
          sel.export_to_csv(mixd, filename)
          
    def export_to_csv(self, mixd, filename):
      '''
      export to csv
      
      '''
        if isinstance(mixd, (pd.DataFrame,)):
            mixd.export_to_csv(filename+'.xlsx')
        else:
            if self.names is None:
                self.names = ['x' + str(n) for n in range(self.nfact)]
            mixd = pd.DataFrame(mixd, columns=self.names)    
        
        
    def export_to_excel(self, mixd, filename):
      '''
      export to excel
      
      '''
        if isinstance(mixd, (pd.DataFrame,)):
            mixd.export_to_excel(filename+'.xlsx')
        else:
            if self.names is None:
                self.names = ['x' + str(n) for n in range(self.nfact)]
            mixd = pd.DataFrame(mixd, columns=self.names)    
    
    def d_optimal(self):
        """
        """
        ...

'''
    def plot_experiments(self):
        """
        """
        if self.with_mixd:
            fig = ff.create_ternary_contour(x, y, pole_labels=self.names, interp_mode='cartesian', colorscale='Viridis', showscale=True, ncontours=ncontours)
            fig.show()
'''

