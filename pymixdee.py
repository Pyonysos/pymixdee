"""
The pymixdee package is designed to help the scientist, statistician, etc., to construct appropriate mixture designs.
"""

import numpy as np
import pandas as pd
from typing import Union
from sympy.utilities.iterables import multiset_permutations

import matplotlib.pyplot as plt


class MixD:
    """
    MixD class:
        attributes:
            nfact:              int, number of factors
            names:              str, list. names of the factors 
            row_limit:          int, maximum number of experiments
        methods:
            dirichlet:          
            centroid_simplex:
            lattice_simplex:
            optimal_mixd:
            add_lower_constraints:
            shuffle:
    """


    def __init__(self, nfact : int=None, factor_names : list = None):
        if (nfact is None) & (factor_names is None):
            raise ValueError('Missing arguments. Expecting value for either nfact or factor_names')

        self.nfact = nfact if nfact is not None else len(factor_names)
        self.names = factor_names if factor_names is not None else ['x'+str(n) for n in range(len(self.nfact))]
        self.with_df_format = True
        
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
    
    def __permutations(self, mixd):
        '''
        permutates mixture proportions to cover evenly the experimental domain
        '''
        for row in mixd:
            permutations = np.array(list(multiset_permutations(row)))
            mixd = np.concatenate((mixd, permutations), axis = 0)
        return mixd

    def __remove_duplicates(self, mixd):
        duplicates = []
        for i in range(mixd.shape[0]):
            for j in range(i+1, mixd.shape[0]):
                if np.all(np.isclose(mixd[i], mixd[j])):
                    duplicates.append(j)
        return mixd[[n for n in range(mixd.shape[0]) if n not in duplicates]]
        
    def df_format(f):
        '''
        decorator converting numpy array to dataframe
        '''
        def inner(self, *args, **kwargs):
            mixd_array = f(self, *args, **kwargs)
            if self.with_df_format == False:
                return mixd_array
            elif self.names:
                return pd.DataFrame(mixd_array, columns=self.names)
            else:
                return mixd_array
        return inner

    @df_format
    def shuffle(self, mixd):
        '''
        method to shuffle the experiments within mix designs
        '''
        mixd = np.array(mixd)
        np.random.default_rng().shuffle(mixd)
        return mixd

    '''
    =================================================================================================
                                    CLASSICAL MIXTURE DESIGNS
    =================================================================================================
    '''   

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
    
    @df_format
    def simplex_centroid(self, ndegree=2, ncenter=1, lower: list= None):
        """
        generate centroid simplex
        ndegree:    int, default 2, number of degree of design
        ncenter:    int, default 1, number of central points to be added
        lower:      list of float, default None, lower constraints to each factor
        """
        ndegree = self.nfact - 1 if ndegree >= self.nfact else ndegree

        base = np.tri(ndegree, self.nfact)
        base /= np.sum(base, axis=1).reshape(-1,1)
        
        base = self.__permutations(base)
        base = self.__remove_duplicates(base)

        if ncenter:
            mixd = self.__add_center_points(base, ncenter)

        print('number of experiments =' , mixd.shape[0])
    
        if lower is not None:
            mixd = self.add_lower_constraints(mixd, lower)
            
        return mixd
        
    @df_format
    def simplex_lattice(self, ndegree=2, ncenter=1, lower: list=None):
        """
        generate Scheffe design
        ndegree:    int, default 2, number of degree of design
        ncenter:    int, default 1, number of central points to be added
        lower:      list of float, default None, lower constraints to each
        """
        
        lattice = np.linspace(0, 1, ndegree, endpoint=False).reshape(-1,1)
        base = np.concatenate((lattice, 1-lattice, np.zeros((lattice.shape[0], self.nfact-2))), axis=1)
        
        #generate multiset_permutations
        base = self.__permutations(base)
        
        #removing duplicates
        base = self.__remove_duplicates(base)

        if ncenter > 0:
            mixd = self.__add_center_points(base, ncenter)
        print('number of experiments =' , mixd.shape[0])

        if lower is not None:
            mixd = self.add_lower_constraints(mixd, lower)
        return mixd


    '''
    =================================================================================================
                                            OPTIMAL DESIGNS
    =================================================================================================
    '''   
    
    def __determinant_infmat(self, mat: np.array):
        '''
        compute determinant of the information matrix:  |X'X|
        '''
        return np.linalg.det( self.__information_matrix(mat) )
    
    def __information_matrix(self, mat: np.array):
        return (mat.T).dot(mat)

    def __dispersion_matrix(self, mat):
        return np.linalg.inv( self.__information_matrix(mat) )
    
    def _d_efficiency(self, mat: np.array):
        '''
        minimiser determinant de la matrice de dispersion (X'X)^-(-1) <=> maximise la matrice d'information (X'X)
        d_eff = 100*det(X'X)^(1/nfact)/ntrial
        '''

        det = self.__determinant_infmat(mat)
       
        return 100 * det**(1/mat.shape[1]) / mat.shape[0]


    def _a_efficiency(self, mat: np.array):
        '''
        minimiser la moyenne de la variance des coefficients de la matrice de dispersion (X'X)^-(-1)
        a-efficiency = (100*nfact)/trace[ ntrial*(X'X)^(-1)]
        '''
        a = mat.shape[1] * self.__dispersion_matrix(mat)
        
        return 100*mat.shape[0] / np.trace(a)
    
    def _g_efficiency(self):
        '''
        minimiser la variance de pr??diction en trouvant les exp??riences qui pr??voient avec le plus de pr??diction
        g_eff = 100 * nfact / (ntrial*d)
        d= maximmum value of x( X^(-1).X )x' (where x is a row of the mixd)
        '''
        ...
    
    @df_format
    def fedorov_algorithm(self, ntrial:int, criterion:str='d', ntest:int=50):
        '''
        1. d??finir un grand nombre d'exp??riences N
        2. d??finir le nombre d'essais n ?? r??aliser
        3. tirer au hasard n exp??riences parmi les N
        4. calculer le crit??re d'optimisation (det(X'X) par exemple)
        5. sortir au hasard une des n exp??riences et ajouter une des N-n
        6. recalculer le crit??re -> si il augmente conserver cet ??change sinon annuler l'echange
        7. it??rer jusque convergence
        '''
        print(f'criterion: {criterion} efficiency')
        #N = np.array(self.dirichlet(ntest*10))
        N = np.array(self.simplex_centroid(ndegree=3, ncenter=1))
        n = N[:ntrial, :]
        convergence_history = []
        for _ in range(ntest):
            i = np.random.default_rng().choice(range(ntrial))
            j = np.random.default_rng().choice(range(ntrial, N.shape[0]))
        
            m = n.copy()
            m[i] = N[j]

            f = getattr(self, f'_{criterion}_efficiency')
            convergence_history.append(self._d_efficiency(n))
            
            if f(m) > f(n):
                n = m

            #n, det = eval(f'self.__{criterion}_efficiency', (n, m))

        plt.plot(range(len(convergence_history)), convergence_history)
        return n

    def genetic_algorithm(self, ntrial : int, criterion: str='d', ntest: int=50):
        '''
        selection
        croisement
        mutation
        '''
        
        N = np.array()



    def optimal_mixd(self, designs: tuple=None, ntrial: int= 20, criteria: str='d', algorithm: str='fedorov', *args, **kwargs):
        """
        plan tq V(b) = (X'X)^(-1) * sig?? /// (X'X)^(-1) d??pend du plan vs. sig?? d??pend des r??sultats 

        => minimiser (X'X)^(-1) 
        Algorithme de fedorov
        """
        ...

    '''
    =================================================================================================
                                            EXPORT METHODS
    =================================================================================================
    '''
    def export(self, mixd:pd.DataFrame, filename:str, extension: str='xlsx'):
        if extension in ('xlsx', 'excel'):
          self.export_to_excel(mixd, filename)
        elif extension == 'csv':
          sel.export_to_csv(mixd, filename)
          
    def export_to_csv(self, mixd: np.array, filename: str):
        '''
        export to csv
        '''
        if isinstance(mixd, (pd.DataFrame,)):
            mixd.export_to_csv(filename+'.csv')
        else:
            if self.names is None:
                self.names = ['x' + str(n) for n in range(self.nfact)]
            mixd = pd.DataFrame(mixd, columns=self.names)
            self.export_to_csv(mixd, filename)
        
        
    def export_to_excel(self, mixd: np.array, filename: str):
        '''
        export to excel
        '''
        if isinstance(mixd, (pd.DataFrame,)):
            mixd.export_to_excel(filename+'.xlsx')
        else:
            if self.names is None:
                self.names = ['x' + str(n) for n in range(self.nfact)]
            mixd = pd.DataFrame(mixd, columns=self.names)
            self.export_to_excel(mixd, filename)