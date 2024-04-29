# -*- coding: utf-8 -*-
"""
Container class for root searching capability of pc_saft interface
contain definition of flattened function for 1,2,3 density root of P(density,temperature) function
"""
__author__ = 'David Celný'
__version__ = '0.0.1'

# Built-in/Generic Imports
# […]
from typing import Optional, Callable

# Libs
from scipy.optimize import newton

# Own modules

class d_root:
    is_silent = False # turn off printing in get_root functions

    def __init__(self, p_function:Callable, target_pressure:float, precision_digits:int=8):
        self.get_p = p_function # has to be the pc_saft_interface type of function returning: p,pp,ppp =f(rho)
        self.target_pressure = target_pressure
        # to bypass repeated calls of function last calculated values are remebered and reused if possible
        self.rho = None
        self.p = None # pressure, output from the eos 
        self.pp = None # denisty derivative of pressure, output from the eos 
        self.ppp = None # second density derivative of pressure, output from the eos 
        # for measuring of number of eos calls
        self.pcall = [0,0,0]
        self.prec_digits = precision_digits
        self.maxiter = 50

    def _clean_root(self, root:float, digits:Optional[int]=None) -> float:
        """ makes sure to get rid of the .99999 noise in numbers """
        if digits is None:
            return round(root, ndigits=self.prec_digits)
        else:
            return round(root, ndigits=digits)

    def get_root0(self, est:float)-> Optional[float]:
        """function that returns the first density root of the p_function"""

        def f_root0(x):            
            self.rho = x
            self.p, self.pp, self.ppp = self.get_p(density=x)
            self.pcall[0] +=1            
            return self.p - self.target_pressure

        def df_root0(x):            
            if x == self.rho:
                pass
            else:
                self.rho = x
                self.p, self.pp, self.ppp = self.get_p(density=x)
                self.pcall[0] +=1
                #NOTE this sets the self.d_iter
            return self.pp
        
        def df2_root0(x):
            if x == self.rho:
                pass
            else:
                self.rho = x
                self.p, self.pp, self.ppp = self.get_p(density=x)
                self.pcall[0] +=1
                #NOTE this sets the self.d_iter
            return self.ppp

        try:
            d_root0, optim_report = newton(f_root0, est,
                                           fprime=df_root0,
                                           fprime2=df2_root0,
                                           full_output=True,
                                           maxiter = self.maxiter,
                                           tol=10**(-self.prec_digits-1))
            
            if optim_report.converged is False:            
                print(f"Unable to converge: {optim_report}")
                d_root0 = None
            else:
                # print(d_root0)
                # print(optim_report)
                d_root0 = self._clean_root(d_root0)
                # print(f"root0: {d_root0}, {self.pcall[0]}")                
            return d_root0

        except RuntimeError as e:
            if self.is_silent is False:
                print("Physical warning: you may be touching un/meta stable conditions where the eos does not work!")
                print(e)
            return None

    def get_root1(self, est:float, d_root0:Optional[float])-> Optional[float]:
        """function that returns the second density root of the p_function"""

        if d_root0 is None:
            if self.is_silent is False:
                print("Not possible to search for more roots: conditions are likely over-ctitical")
            return None
        
        def f_root1(x, y0):            
            self.rho = x
            self.p, self.pp, self.ppp = self.get_p(density=x)
            self.pcall[1] +=1
            return (self.p - self.target_pressure)/(x - y0)

        def df_root1(x, y0):            
            if x == self.rho:
                pass
            else:
                self.rho = x
                self.p, self.pp, self.ppp = self.get_p(density=x)
                self.pcall[1] +=1
            # NOTE in both cases self.rho == x, so it is used same as the self.p, self.pp, self.ppp
            # NOTE analytic derivative of collapsed(divided by first root) function
            # WOLFRAM ALPHA: d/dx((f(x) - p)/(x - y0)) = ((x - y0) f'(x) - f(x) + p)/(x - y0)^2
            return ((self.rho - y0)*self.pp - self.p + self.target_pressure)/(self.rho - y0)**2
        
        def df2_root1(x, y0):
            if x == self.rho:
                pass
            else:
                self.rho = x
                self.p, self.pp, self.ppp = self.get_p(density=x)
                self.pcall[1] +=1
            # NOTE in both cases self.rho == x, so it is used same as the self.p, self.pp, self.ppp
            # NOTE analytic derivative of collapsed(divided by first root) function
            # WOLFRAM ALPHA: d^2/dxdx((f(x) - p)/(x - y0)) = ((x - y0)^2 f''(x) + 2 (y0 - x) f'(x) + 2 f(x) - 2 p)/(x - y0)^3            
            return ((self.rho - y0)**2 *self.ppp + 2*(y0 - self.rho)*self.pp + 2*self.p - 2*self.target_pressure)/(self.rho - y0)**3
        try:
            d_root1, optim_report = newton(f_root1, est,
                                            args=(d_root0,),
                                            # fprime=df_root1,
                                            # fprime2=df2_root1,
                                            full_output=True,
                                            maxiter = self.maxiter,
                                            tol=10**(-self.prec_digits-1))
            
            if optim_report.converged is False:            
                # NOTE - can be case of only single root (i.e. vapor or liquid)
                d_root1 = None
            else:                
                d_root1 = self._clean_root(d_root1)
                # print(f"root1: {d_root1}, {self.pcall[1]}") # DEBUG
            return d_root1

        except RuntimeError as e:
            # NOTE - can be case of only single root (i.e. vapor or liquid)
            if self.is_silent is False:
                print("Physical warning: you may be touching unstable conditions where the eos does not work!")
                print(e)
            return None

    def get_root2(self, est:float, d_root0:Optional[float], d_root1:Optional[float])-> Optional[float]:
        """function that returns the third density root of the p_function"""

        if d_root0 is None:
            if self.is_silent is False:
                print("Not possible to search for more roots: conditions are likely over-ctitical")
            return None
        
        if d_root1 is None:
            if self.is_silent is False:
                print("Not possible to search for more roots: conditions are likely unstable")
            return None

        def f_root2(x, y0, y1):            
            self.rho = x
            self.p, self.pp, self.ppp = self.get_p(density=x)
            self.pcall[2] +=1
            return (self.p - self.target_pressure)/((x - y0)*(x - y1))

        def df_root2(x, y0, y1):            
            if x == self.rho:
                pass
            else:
                self.rho = x
                self.p, self.pp, self.ppp = self.get_p(density=x)
                self.pcall[2] +=1
            # NOTE in both cases self.rho == x, so it is used same as the self.p, self.pp, self.ppp
            # NOTE analytic derivative of collapsed(divided by first root) function            
            # WOLFRAM ALPHA: d/dx((f(x) - p)/((x - y0) (x - y1))) = ((x - y0) (x - y1) f'(x) + f(x) (-2 x + y0 + y1) + p (2 x - y0 - y1))/((x - y0)^2 (x - y1)^2)
            return ((self.rho - y0)*(self.rho - y1)*self.pp + self.p*(-2*self.rho + y0 + y1) + self.target_pressure*(2*self.rho - y0 - y1))/((self.rho - y0)**2 *(self.rho - y1)**2)
        
        def df2_root2(x, y0, y1):
            if x == self.rho:
                pass
            else:
                self.rho = x
                self.p, self.pp, self.ppp = self.get_p(density=x)
                self.pcall[2] +=1
            # NOTE in both cases self.rho == x, so it is used same as the self.p, self.pp, self.ppp
            # NOTE analytic derivative of collapsed(divided by first root) function
            # WOLFRAM ALPHA:d^2/(dx dx)(f(x) - p)/((x - y0) (x - y1)) = -(-(x - y0)^2 (x - y1)^2 f''(x) + 2 (x - y0) (x - y1) f'(x) (2 x - y0 - y1) - 2 f(x) (3 x^2 - 3 x (y0 + y1) + y0^2 + y0 y1 + y1^2) + 2 p (3 x^2 - 3 x (y0 + y1) + y0^2 + y0 y1 + y1^2))/((x - y0)^3 (x - y1)^3)
            return -(-(self.rho - y0)**2 *(self.rho - y1)**2 *self.ppp + 2*(self.rho - y0)*(self.rho - y1)*self.pp*(2*self.rho - y0 - y1) - 2*self.p*(3*self.rho**2 - 3*self.rho*(y0 + y1) + y0**2 + y0*y1 + y1**2) + 2*self.target_pressure*(3*self.rho**2 - 3*self.rho*(y0 + y1) + y0**2 + y0*y1 + y1**2))/((self.rho - y0)**3 *(self.rho - y1)**3)
        
        try:
            d_root2, optim_report = newton(f_root2, est,
                                            args=(d_root0, d_root1),
                                            # fprime=df_root2,
                                            # fprime2=df2_root2,
                                            full_output=True,
                                            maxiter = self.maxiter,
                                            tol=10**(-self.prec_digits-1))
        
            if optim_report.converged is False:            
                # NOTE - can be case of only single root (i.e. vapor or liquid)
                d_root2 = None
            else:
                d_root2 = self._clean_root(d_root2)
                # print(f"root2: {d_root2}, {self.pcall[2]}") # DEBUG
            return d_root2

        except RuntimeError as e:
            # NOTE - can be case of two roots in metastable region 
            if self.is_silent is False:
                print("Physical warning: you may be touching metastable conditions where the eos does not work!")
                print(e)
            return None
        
        