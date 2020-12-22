import thermo
import scipy
import numpy as np
import matplotlib.pyplot as plt

def diff(y,dx):
    out = np.ndarray(np.size(y))

    out[0] = (-3*y[0]+4*y[1]-1*y[2])*(2*dx)
    out[-1] = (3*y[-1]-4*y[-2]+1*y[-3])*(2*dx)
    for i in range(1,len(out)-1):
        out[i] = (y[i+1]-y[i-1])/(2*dx)
    return out

def integ(y,dx,l1_i,end=-1,start=0):
    out = np.zeros(len(y[start:end]))
    a_end = end-start
    a_l1_i = l1_i-start
    out[a_l1_i] = 0
    if l1_i != end:
        for i in range(a_l1_i+1,len(out)):
            out[i] = out[i-1]+y[i-1]*dx+(y[i]-y[i-1])*0.5*dx
    if l1_i != start:
        for i in range(a_l1_i-1,-1,-1):
            out[i] = out[i+1]-y[i+1]*dx-(y[i]-y[i+1])*0.5*dx
    
    return out

def refine(geometry,factor):
    f = scipy.interpolate.interp1d(geometry[:,0], geometry[:,1])
    dx = geometry[1,0]-geometry[0,0]
    x = np.arange(geometry[0,0],geometry[-1,0],dx/factor)
    return np.transpose(np.array([x,f(x)]))


class x_bar:
    #class for effective controur distance
    def __init__(self,R_x,x_loc_ind,dx,x_throat_ind):
        self.r_i = R_x[x_loc_ind]
        self.dx_dy = diff(R_x,dx)
        self.R_x = R_x
        self.dx = dx
        self.x_ti = x_throat_ind

    def compute(self,x_loc_ind,end=-1,start=0):
        a = 1/((x_loc_ind-self.x_ti)*self.dx)
        b = 0.75-a*x_loc_ind*self.dx
        x = np.arange(0,len(self.R_x)*self.dx,self.dx)
        return integ(self.R_x/self.r_i*np.sqrt(self.dx_dy*self.dx_dy+1)*(a*x+b),self.dx,x_loc_ind,end=end,start=start)
        
        

class film:
    def __init__(self):
        #def __init__(self, massflow, V_channel,V_chamber, x_loc_ind, inlet_angle,R,P_channel,T_channel,P_chamber,fluid,mixture,wt_channel):
        """
        self.fluid = thermo.Mixture(fluid, ws=mixture,  T=T_channel[x_loc_ind], P=P_channel[x_loc_ind])
        self.wt = wt_channel[x_loc_ind]
        self.massflow = massflow
        self.pressuredrop = (P_channel[x_loc_ind] - 1/2*self.fluid.rho*V_channel[x_loc_ind]**2)-P_chamber
        self.inletangle = inlet_angle
        self.V_chamber = V_chamber[x_loc_ind]
        self.R = R[x_loc_ind]
        """
        self.fluid = thermo.Mixture(['c2h5oh','h2o'], [0.9,0.1],  T=350, P=65e5)
        self.wt = 0.6e-3
        self.massflow = 0.2
        self.pressuredrop = 15e5
        self.inletangle = 30/180*np.pi
        self.V_chamber = 200

        self.R = 25e-3

    

    def xi1c(self, Re):
        # curfve fit from "Liquid Rocket Trhust Chambers" page 48 realting inlet efficiency to Reynolds number
        Re = np.log10(Re)
        return 3.55378*np.exp(-Re*0.647016) - 0.103358

    def xiinlet(self):
        return 0.5 + 1.2/np.pi*self.inletangle         # min 0.5 for coaxial flow before injection, max 0.9 for flow at pi/6 rad realtive to faceplate 
    
    def injector(self, maxiter=100, tol=1e-6):


        xiinlet = self.xiinlet()
        self.xi = xiinlet
        
        def func(diameter):
            self.length = self.wt/np.sin(self.inletangle)
            mdt_hole = self.massflow* diameter/(self.R*2*np.pi)
            self.mu = 1/np.sqrt(1 + self.xi)
            diameter = 0.95*mdt_hole**0.5 * self.mu**(-0.5) * (self.fluid.rho*self.pressuredrop)**(-0.25)
            self.vel = 1.273*mdt_hole*self.fluid.rho**(-1)*diameter**(-2)
            self.Re = self.fluid.rho*self.vel*diameter/self.fluid.mu
            lam = 0.3164*self.Re**(-0.25)
            friction = lam*self.length/diameter
            xi1c = self.xi1c(self.Re)

            self.xi = xiinlet + xi1c + friction 
            newmu = 1/np.sqrt(1 + self.xi)

            return abs(self.mu - newmu)

        sol = scipy.optimize.root_scalar(func, x0 = 0.1e-3,x1 = 0.2e-3, method='secant')
        self.diameter = sol.root
        #print(self.diameter)

if __name__ == "__main__":
	geometry = np.genfromtxt('sparrow_50bar.txt', delimiter='', dtype=None, skip_header = 13) / 1000
	rg = refine(geometry,10)
	rx = rg[:,1]
	dx = rg[1,0]-rg[0,0]
	x = x_bar(rx,350,dx,500)

	plt.plot(x.compute(350,700,340))
	plt.show()
