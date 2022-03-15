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


class DiffNotation:
    '''Creates a DifNotation object

    Parameters
    ----------
    equation : sympy.core.add.Add
       A sympy made equation.
    *args : sympy.core.symbol.Symbol
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