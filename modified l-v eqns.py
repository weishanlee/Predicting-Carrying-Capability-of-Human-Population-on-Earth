# -*- coding: utf-8 -*-
"""
Corresponding author: weishan_lee@yahoo.com

Description
---------------------------
Earth's Carrying Capability of Earth Human Population. 
This code provides details published in Ref[1].
References:
[1]    
[2] Rein Taagepera,"A world population growth model: Interaction with Earth's 
    carrying capability and technology in limited space", Technological 
    Forecasting & Social Change 82 (2014) 34-41.
[3] https://www.kaggle.com/theworldbank/global-population-estimates

"""
#print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

nsteps = 2500    # time step (month)

timeStep = np.arange(0, nsteps, 1)
for i, val in enumerate(timeStep):
    timeStep[i] = 2019 + val/12.0
#%% Modified Lotka-Volterra Equations
# fitting parameters for the t-function
A = 3.83e9
B = 1.28
D = 1980.0
tau = 22.9 # years
M = 0.70

# fitting parameters for the Lotka-Volterra Equation
alpha_ = 23.0457
beta = alpha_ / (np.e)**2 # 3.11889633

gamma_ = 23.0457
delta = gamma_ / (np.e)**2 # 3.11889633
eta = (np.e)**2
zeta = 25 

def EFunc(t):
    value = np.exp( (D-t)/tau )
    return value

def Resources(t):
    value = eta + M * EFunc(t) / ( delta * tau * (B + EFunc(t)) * np.log( B + 
            EFunc(t) ) ) 
    return  value

def x2(t):
    value = eta - M * EFunc(t) / ( beta * tau * (B + EFunc(t)) * np.log( B + 
            EFunc(t) ) ) 
    return  value

def tFunc(t): # t-function
    value = A / ( np.log( (B + EFunc(t)) ) )**M
    return value

def gamma(t):
    numerator =  M * (np.log(B+EFunc(t)))**(M-1)*EFunc(t)\
                 * (np.log(B+EFunc(t)) * B - EFunc(t))
    denominator = tau*(-eta*tau*beta*(B+EFunc(t))*(np.log(B+EFunc(t)))**(M+1)
                      +(np.log(B+EFunc(t)))**M*EFunc(t)*M+A*zeta*
                      (eta*tau*beta*(B+EFunc(t))*np.log(B+EFunc(t))-M*EFunc(t))
                      )*(B+EFunc(t))
    value = numerator/denominator
    return value

def alpha(t):
    numerator =  M * EFunc(t) * (np.log(B+EFunc(t)))**(M-1)\
                 * ( -np.log(B+EFunc(t)) * B + EFunc(t))
    denominator = tau*(-eta*tau*delta*(B+EFunc(t))*(np.log(B+EFunc(t)))**(M+1)
                      -(np.log(B+EFunc(t)))**M*EFunc(t)*M+A*zeta*
                     (eta*tau*delta*(B+EFunc(t))*np.log(B+EFunc(t))+M*EFunc(t))
                      )*(B+EFunc(t))
    value = numerator/denominator
    return value

def KFunc1(t):
    #value = (eta / np.e)**alpha_ * (1/(zeta*np.e))**gamma(t)
    x = tFunc(t)
    y = x2(t)
    delta = zeta * gamma(t)
    value = y**alpha * np.exp(-beta*y) * x**gamma(t) * np.exp(-delta*x)    
    return value

def KFunc2(t):
    #value = (eta / np.e)**gamma_ * (1/(zeta*np.e))**alpha(t)
    x = Resources(t)
    y = tFunc(t)
    beta = zeta * alpha(t)
    delta = gamma_ / (np.e)**2
    value = y**alpha(t) * np.exp(-beta*y) * x**gamma_ * np.exp(-delta*x)    
    return value

# plot T-Function and recorded data of fitting T-Function
plt.figure("T-Function")
ax = plt.gca()
ax.set_xlabel("Time (year)",size = 16)
ax.set_ylabel("Population Size",size = 16)
plt.grid()    

xvalsRef1 = [400,600,800,1000,1100,1200,1300,1400,1500,1600,1700,1750,1800,
             1850,1900,1920,1940,1950,1960,1970,1980,1990,2000,2010]
yvalsRef1 = [198,214,235,281,310,398,396,362,457,544,635,771,941,1242,1639,
             1905,2313,2526,3035,3667,4442,5278,6021,6861] 
for i, val in enumerate(yvalsRef1):
    yvalsRef1[i] = val * 1e6

xvalsTimeStep = np.arange(-2,2500, 1)
yvalsx1 = tFunc(xvalsTimeStep)

plt.plot(xvalsRef1, yvalsRef1,'ro',xvalsTimeStep, yvalsx1,'k-',lw=1)

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5)#number of minor intervals per major inteval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
lines=ax.get_lines()
lines[0].set_label("Data from Table 1 in Ref[2]")
lines[1].set_label("T-function")
ax.legend(loc='center left', shadow=True, fontsize='large')
plt.show()
plt.savefig("./tfunction.png")

# plot all together
plt.figure("Population Size vs. Time (year)")
ax = plt.gca()
ax.set_xlabel("Time (year)",size = 16)
ax.set_ylabel("Population Size",size = 16)
plt.grid()

xvalsRef2 = np.arange(1962,2051,10)

yvalsRef2 = [3127961482,3192794384,3258201476,3324951621,3394864530,3464439525,
3534821115,3609383725,3684765870,3762198347,3838924951,3914857611,	
3991430917,4066267984,4139151082,4211781677,4285609387,4361295248,	
4437690434,4515764583,4596813158,4678525765,4759982360,4843067309,	
4928822143,5016798785,5105701987,5194731380,5284886348,5372078249,
5456141249,5541075501,5624840414,5709757338,5792568347,5875398158,
5957237460,6038067278,6118075293,6197638117,6276824418,6356259574,	
6436346998,6517020798,6598421257,6680423047,6763745673,6847214549,	
6930656699,7012843635,7097400665,7182860115,7268986176,7355220412,	
7442135578,7524453000,7606102000,7686852000,7766687000,7845612000,	
7923611000,8000663000,8076739000,8151786000,8225788000,8298763000,	
8370716000,8441688000,8511707000,8580819000,8649017000,8716283000,	
8782656000,8848138000,8912769000,8976588000,9039596000,9101784000,
9163183000,9223737000,9283410000,9342117000,9399793000,9456429000,	
9511938000,9566349000,9619629000,9671774000,9722319000]

yvalsRef2 = [yvalsRef2[i] for i in range(len(yvalsRef2)) if i%10 == 0]

yvalsKFunc2 = KFunc2(xvalsTimeStep)

plt.plot(xvalsRef1, yvalsRef1,'ro',
         xvalsRef2, yvalsRef2,'mo',
         xvalsTimeStep, yvalsx1,'k-',
         xvalsTimeStep, yvalsKFunc2,'g-.',lw=1)

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5)#number of minor intervals per major inteval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
lines=ax.get_lines()
lines[0].set_label("Data from Table 1 in Ref[2]")
lines[1].set_label("Data from Ref[3]")
lines[2].set_label("T-function")
lines[3].set_label("Carrying Capability K(t)")
ax.legend(loc='center left', shadow=True, fontsize='large')
plt.xlim(1500, 2200)
plt.show()
plt.savefig("./allTogether.png")

# plot Resources
plt.figure("x1(t)")
ax = plt.gca()
ax.set_xlabel("Time (year)",size = 16)
ax.set_ylabel("Natural resources $x_{1}(t)$",size = 16)
plt.grid()

yvals = Resources(xvalsTimeStep)

plt.plot(xvalsTimeStep,yvals,'k-',lw=2)

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5)#number of minor intervals per major inteval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis

plt.show()
plt.savefig("./x1.png")

# plot x2
plt.figure("x2(t)")
ax = plt.gca()
ax.set_xlabel("Time (year)",size = 16)
ax.set_ylabel("Lethal factors $x_{2}(t)$",size = 16)
plt.grid()

yvals = x2(xvalsTimeStep)

plt.plot(xvalsTimeStep,yvals,'k-',lw=2)

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5)#number of minor intervals per major inteval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis

plt.show()
plt.savefig("./x2.png")

# plot alpha(t)
plt.figure("alpha(t)")
ax = plt.gca()
ax.set_xlabel("Time (year)",size = 16)
ax.set_ylabel("$\\alpha(t)$",size = 16)
plt.grid()

yvals = alpha(xvalsTimeStep)

plt.plot(xvalsTimeStep,yvals,'k-',lw=2)

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5)#number of minor intervals per major inteval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis

plt.show()
plt.savefig("./alpha.png")

# plot gamma(t)
plt.figure("gamma(t)")
ax = plt.gca()
ax.set_xlabel("Time (year)",size = 16)
ax.set_ylabel("$\gamma(t)$",size = 16)
plt.grid()

yvals = gamma(xvalsTimeStep)

plt.plot(xvalsTimeStep,yvals,'k-',lw=2)

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5)#number of minor intervals per major inteval
minorLocatorY = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis

plt.show()
plt.savefig("./gamma.png")