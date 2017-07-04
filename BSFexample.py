
# coding: utf-8

# In[ ]:

# Python implementation of the Black-Scholes functions.
# BE AWARE of the fact that if you assign spot and strike values as integers, 
# then the critical "moneyness" quantity spot/strike will be calculated in 
# integer arithmetic.  In this case, if spot < strike, the quotient will be 
# 0 and the logarithm function will quit on you.
#
# Therefore, use decimal points in the inputs that are not integers.
# In Summary: BSF_call(50, 0.01, 0.25, 1.0, 45, 0) will give an error
#             BSF_call(50.0, 0.01, 0.25, 1.0, 45.0, 0) will not
  
#from Normal import *

import math

def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    t = 1. / (1. + 0.5*z)
    r = t * math.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
    	t*(.09678418+t*(-.18628806+t*(.27886807+
    	t*(-1.13520398+t*(1.48851587+t*(-.82215223+
    	t*.17087277)))))))))
    if (x >= 0.):
    	return r
    else:
    	return 2. - r

def Normal(x):
    return 1. - 0.5*erfcc(x/(2**0.5))

def NormalDen(x):
    return (1.0/(2.0*3.141592653589793238)**0.5)*math.exp(-0.5*x*x)

def BSF_call(spot, rate, vol, T, strike, div):
    sqrtT = math.sqrt(T)        # compute this on its own since it is used in more than one place
    d1 = ( math.log(spot/strike) + (rate - div + 0.5*vol*vol)*T )/(vol*sqrtT)
    #sqrtT = sqrt(T)        # compute this on its own since it is used in more than one place
    #d1 = ( log(spot/strike) + (rate - div + 0.5*vol*vol)*T )/(vol*sqrtT)
    d2 = d1 - vol*sqrtT
    return spot*math.exp(-div*T)*Normal(d1) - strike*math.exp(-rate*T)*Normal(d2)
    
def BSF_put(spot, rate, vol, T, strike, div):
    sqrtT = math.sqrt(T)        # compute this on its own since it is used in more than one place
    d1 = ( math.log(spot/strike) + (rate - div + 0.5*vol*vol)*T )/(vol*sqrtT)
    d2 = d1 - vol*sqrtT
    return -spot*math.exp(-div*T)*Normal(-d1) + strike*math.exp(-rate*T)*Normal(-d2)
    
def BSF_Test(spot, rate,vol,T,strike, div, gap):   # check of put-call parity over a range of spot prices
    s = spot - gap
    ds = 2*gap/10.0
    disc = math.exp(rate*T)   # the reciprocal of the discount factor
    while s <=spot + gap:
        SmCpP =  s - BSF_call(s, rate, vol, T, strike, div) + BSF_put(s, rate, vol, T, strike, div)  
        s += ds
        print math.exp(rate*T)*SmCpP   # we should see the strike value each time
        
# Also checked with on-line calculator
        
# A function-object version of BSF in which volatility is the variable.  
class BSF(object):    
    def __init__(self, spot, rate, T, strike, div, isCall):  #the "constructor of the object"
        self.spot = spot   
        self.rate = rate   
        self.T = T   
        self.sqrtT = math.sqrt(self.T)
        self.strike = strike   
        self.div = div   
        self.isCall = isCall    #set this equal to 0 for a put and to 1 for a call
        self.div_disc = math.exp(-div*self.T)

    def __call__(self, vol):    #this is like the C++ "double operator()(double vol)"  
        if(self.isCall): return BSF_call(self.spot, self.rate, vol, self.T, self.strike, self.div)
        else: return BSF_put(self.spot, self.rate, vol, self.T, self.strike, self.div)
        
    def vega(self, vol):  #this is the derivative of the BSF with respect to volatility; needed in Newton's method 
         d1 = ( math.log(self.spot/self.strike) + (self.rate - self.div + 0.5*vol*vol)*self.T )/(vol*self.sqrtT)
         return self.spot*self.div_disc*NormalDen(d1)*self.sqrtT

    def impVol(self, mktP, v_tmp=0.1, tol=1.0e-6):
        iter_num = 0
        iter_max = 50
        while abs(mktP - self.__call__(v_tmp)) >= tol and iter_num < iter_max:
            v_new = v_tmp - (self.__call__(v_tmp) - mktP)/self.vega(v_tmp)
            v_tmp = v_new
            iter_num = iter_num +1
        else:
            if iter_num >= iter_max:
               print 'Max iteration reached before finding implied vol'
               print 'Pricing error is: '+ str(abs(mktP - self.__call__(v_tmp)))
            print 'Iteration number is: '+str(iter_num)
            return v_tmp

# the BSF class definition is finished
        
def price(spot, rate, vol, T, strike, div, isCall):  #a round about way to evaluate the BSF    
    f = BSF(spot, rate, T, strike, div, isCall)
    return f(vol)

def volTest(spot, rate, vol0, vol1, T, strike, div, isCall, num_steps): #will evaluate BSF with vol values between vol0 and vol1, taking num_steps
    f = BSF(spot, rate, T, strike, div, isCall)
    s = vol0
    step = (vol1 - vol0)/num_steps
    while s <= vol1:
        print f(s)
        s = s + step

if __name__ == "__main__":
    # Set parameters for options
    S_0 = 15247.92
    rf_rate = 0.025
    mature_T = 32.0/247.0
    strike = 17000.0
    q_div = 0
    isCall = 1  # value 1: call option;  value 0: put option
    
    ## instantiate a call option: call_opt1
    call_opt1 = BSF(S_0, rf_rate, mature_T, strike, q_div, isCall)
    # Set the volatility of the stock price distribution model (log-normal model)
    vol_sigma = 0.22
    print 'Option price is: '+str(call_opt1.__call__(vol_sigma))
    print 'Option vega is: '+str(call_opt1.vega(vol_sigma))
    mkt_P = 54.006
    imp_vol = call_opt1.impVol(mkt_P,0.15)
    print 'Implied vol is: '+ str(imp_vol)
    


