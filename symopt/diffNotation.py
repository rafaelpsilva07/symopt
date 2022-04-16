'''
Differentiation Notation
========================

This module  computes differentiation notation for functions of several variables
Based on page 90 of reference 1.

References
==========
[1] ARORA, Jasbir S. Introduction to Optimum Design. Elservie. 2nd ed.

'''

import sympy as sp

#from symopt.optimization_functions.goldenSearch import goldenSearch


import sys
path = 'C:/Users/Pichau/repositories_C/symopt/symopt'
sys.path.append(path)
from optimization_functions.goldenSearch import goldenSearch


class DiffNotation:
    '''Creates a DifNotation object

    Parameters
    ----------
    equation : sympy.core.add.Add
       A sympy made equation.
    args : sympy.core.symbol.Symbol
       Symbolic variable that will be differentiated.

    Examples
    --------
    >>> import sympy as sp
    >>> x1, x2 = sp.symbols(['x1', 'x2'])
    >>> fx1x2 = x1**2 + x1*x2 + x2**2
    >>> DiffNotation(fx1x2, x1, x2)
    '''

    def __init__(self, equation, *args):
        self.equation = equation
        self.args = args
        self.xi_pk = None # parametrization iteration k
        self.grad_k = None # gradient of iteration k

    def gradient(self):
        """ Get the gradient of the equation

        Parameters
        ----------
        self : self
            Instance of DifNotation
        
        Returns
        -------
        result : sp.Matrix
            The gradient of the equation in vector notation.
        """
        
        result = [
            sp.diff(self.equation, variable) for variable in self.args
        ]
        return sp.Matrix(result)
        
    def hessian(self):
        """ Get the gradient of the equation

        Parameters
        ----------
        self : self
            Instance of DifNotation
        
        Returns
        -------
        result : sp.Matrix
            The hessian of the equation in vector notation.
        """

        result = [[sp.diff(self.equation, vari, varj) 
                   for vari in self.args] for varj in self.args]             
        return sp.Matrix(result)

    def _parametrization(self, X0):
        '''Parametrize equation in function of alpha'''

        a = sp.symbols('alpha')
        
        xi_0 = {} # X0 dictionary to apply .subs
        for i, value in enumerate(self.args):
            xi_0[value] = X0[i]

        self.grad_k = -self.gradient().subs(xi_0)

        xi_pk = {} # dictionary to apply .subs
        dk = {}  # directions
        for i, value in enumerate(self.args):
            xi_pk[value] = X0[i] + a * self.grad_k[i]

        fa = self.equation.subs(xi_pk)
        self.xi_pk = xi_pk # creates atribute of parametrization
        self.dk = -self.grad_k # initiate attribute of direction
        return fa

    def optimize(
            self,
            X0, 
            epsilon=0.0001, 
            method='goldenSearch', 
            direction='conjugate', 
            I=2e-4, 
            delta=0.1,
            debug_CG=False,
            debug_goldenSearch=False):
        '''
        Parameters
        ----------
        X0 : iterable
            An iterable containing the initial value
        epsilon : float , default 0.0001
            Error based in gradient norm
        method : string, default 'goldenSearch '
            GoldenSearch is the only method available
        direction : string, default 'conjugate'
            Gradient Conjugate is the only method implemented.
        I : float, default : 2e-4
            Error to apply to the goldenSearch method
        delta : float, default 0.1
            Step to apply to the goldenSearch method
        
        Returns
        -------
        res : dict
            Dict containing the results of optimization.
        '''
        
        a = sp.symbols('alpha')
        error = 100 # large error value to start the while loop
        fa = self._parametrization(X0) # first parametrization
        alpha = goldenSearch(fa, delta, debug_goldenSearch, I)

        k = 0 # conjugate gradient number of directions
        res = {} # dictionary that returns the solution

        for variable in self.args:
            res[variable] = []

        while abs(error) > epsilon:
            xk_n = {} # xk of current iteration
            for i, value in enumerate(self.args):
                xk_n[value] = self.xi_pk[value].subs(a, alpha)
        
            grad_k1 = -self.gradient().subs(xk_n)
            beta_kn = (grad_k1.norm() / self.grad_k.norm())**2
            self.grad_k = grad_k1 # updates self.grad_k

            
            # update parametrization variables and error
            e = 0
            for i, value in enumerate(self.args):
                d_kn = grad_k1[i] + beta_kn * self.dk[i]
                e += d_kn**2
                self.xi_pk[value] = xk_n[value] + a * d_kn
                self.dk[i] = d_kn # updates dk

            fa = self.equation.subs(self.xi_pk) # update parametrization
            alpha = goldenSearch(fa, delta, debug_goldenSearch, I) # update alpha
            error = e ** 0.5 # norm of d_kn
            k += 1

            if debug_CG:
                for variable, value in self.xi_pk.items():
                    print(f'{variable} =  {xk_n[variable]}')
                    res[variable].append(xk_n[variable])

            if k > 500:
                print(f'alpha = {alpha}')
                print(f'xi_pk = {self.xi_pk}')
                print(f'f* = {fa.subs(a, alpha)}')
                return ValueError('Convergence not reached until 500 direction searchs')

        # solution

        res.update({'alpha': alpha,
               'f': fa.subs(a, alpha),
               'k_direction_searchs': k})

        if not debug_CG:
            for i, value in enumerate(self.args):
                res[value] = self.xi_pk[value].subs(a, alpha)
        return res

