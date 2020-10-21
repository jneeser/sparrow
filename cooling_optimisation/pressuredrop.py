import numpy as np
import geom_class as ge
import scipy.optimize
from scipy import interpolate
import time

# TODO un-hardcoade TBC heat flux params
mdt = 2.21 #kg/s
k = 6e-6 #guess
ki = 0.167

def f(D,Rev,implicit=True):
    #colebrook equation
    if implicit:
        def func(f):
            return -2*np.log10(k/(3.7*D)+2.51/(Rev*np.sqrt(f)))-1/np.sqrt(f)
        res = scipy.optimize.root_scalar(func,bracket=[1e-10,1]).root
    else:
    #explisict aproximation from https://www.sciencedirect.com/science/article/abs/pii/S0920410513003495
        G = 6.0173 / (Rev*(0.07*k/D+Rev**-0.085))+k/(D*3.71)
        res = ((2.51/Rev+1.1513*G) / (G - k/(D*3.71)-2.3026*G*np.log10(G)))**2
    return res

def Kt(Dd,rd):
    #interpolation from https://web.njit.edu/~sengupta/met%20301/Stress%20Concentration.pdf trandeline in excel from 
    interp_const = np.array([[1.01,1.02,1.03,1.05,1.07,1.1,1.2,1.3,2,3,6],
                             [0.9327,0.9657,0.9896,0.9995,1.0149,1.0093,0.9884,0.9957,0.9811,0.9354,0.9141],
                             [-0.169,-0.179,-0.188,-0.199,-0.203,-0.217,-0.24,-0.256,-0.281,-0.318,-0.346]])

    A = interpolate.interp1d(interp_const[0],interp_const[1],fill_value='extrapolate')
    b = interpolate.interp1d(interp_const[0],interp_const[2],fill_value='extrapolate')
    return A(Dd)*rd**b(Dd)

def deltaP(D,fv,rho,v):
    #darcy wisebach pressure drop per unit length
    return fv*rho*v**2/(D*2)

def Re(D,rho,v, mu):
    #reynold number
    return rho*v*D/mu

class parameters:
    def __init__(self,ri,t,wt1,wt2,rf1,rf2,N):
        self.ri = ri
        self.t = t
        self.ro = ri + self.t + wt1
        self.wt1 = wt1
        self.wt2 = wt2
        self.rf1i = rf1
        self.rf1o = rf1
        self.rf2 = rf2
        self.N = N
        self.wt3 = None

    def update(self,t,wt1,rf1i,rf1o,rf2):
        self.ro = self.ri + t + wt1
        self.t = t
        self.wt1 = wt1
        self.rf1i = rf1i
        self.rf1o = rf1o
        self.rf2 = rf2

class metal:
    def __init__(self,E,k,v,alpha,sig_yield):
        self.E = E
        self.k = k
        self.v = v
        self.alpha = alpha
        self.sig_yield_table = sig_yield
        self.sig_yield = interpolate.interp1d(sig_yield[1],sig_yield[0],fill_value='extrapolate')

class sim:
    def __init__(self, max_wall_temp, T_wall,T_cool,h_cha,q_rad,p_cool,p_cha,y_coordinate,SF=1.5):
        self.max_wall_temp = max_wall_temp
        self.T_wall = T_wall
        self.T_cool = T_cool
        self.h_cha = h_cha
        self.q_rad = q_rad
        self.p_cool = p_cool
        self.p_cha = p_cha
        self.SF = SF
        self.pd_con = 200000/(y_coordinate)


class physics:
    def __init__(self,x,metal,sim, fluid):
        self.par = x
        self.sim = sim
        self.met = metal
        self.record_cross = []
        self.record_geom = []
        self.fluid = fluid

        def geom_update(input):
            self.par.update(input[0],input[1],input[2],input[3],input[4])
            self.N = x.N
            self.channeli , self.Do = geomi(self.par)
            self.channelo , self.Di = geomo(self.par)
            self.dhi = 4*self.channeli.A/self.channeli.circ
            self.dho = 4*self.channelo.A/self.channelo.circ

            Lo =  np.linalg.norm(self.channelo.pointclasslist[0][0]-self.channelo.pointclasslist[1][0])

        def flow_update():
            V_total = mdt/(self.fluid.rho*self.N)
            def func(V_rat):
                self.Vi = V_total*V_rat
                self.Vo = V_total-self.Vi
                self.vi = self.Vi/self.channeli.A
                self.vo = self.Vo/self.channelo.A
                self.Rei = Re(self.dhi,self.fluid.rho,self.vi, self.fluid.mu)
                self.Reo = Re(self.dho,self.fluid.rho,self.vo, self.fluid.mu)
                self.fi = f(self.dhi,self.Rei)
                self.fo = f(self.dho,self.Reo)
                return deltaP(self.dhi,self.fi,self.fluid.rho,self.vi)-deltaP(self.dho,self.fo,self.fluid.rho,self.vo)
            
            scipy.optimize.root_scalar(func,bracket=[0.1,0.9])
            #print(self.Vi/self.Vo)
        
        def stress_constraint(input):
            #print(input)
            geom_update(input)
            flow_update()
            self.t_tbc = 0.15e-3
            self.k_tbc = 2
            self.Pr = self.fluid.Pr
            self.hi = ki/self.dhi*0.023*self.Rei**0.8*self.Pr**0.4
            self.ho =  ki/self.dhi*0.023*self.Reo**0.8*self.Pr**0.4
            self.q = (self.sim.T_wall - self.sim.T_cool + self.sim.q_rad/self.sim.h_cha) / (1/self.sim.h_cha + self.par.wt1/self.met.k + 1/self.hi + self.t_tbc/self.k_tbc)
            Li =  np.linalg.norm(self.channeli.pointclasslist[2][0]-self.channeli.pointclasslist[3][0])

            temp_sigma_steady1 = self.met.E * self.met.alpha * self.q * self.par.wt1 / (2*(1-self.met.v)*self.met.k) - (self.sim.p_cool-self.sim.p_cha)*Li*Li*3/(4*self.par.wt1**2)
            temp_sigma_steady2 = self.met.E * self.met.alpha * self.q * self.par.wt1 / (2*(1-self.met.v)*self.met.k) - (self.sim.p_cool-self.sim.p_cha)*Li*Li/(4*self.par.wt1**2)
            temp_sigma_start = (self.sim.p_cool)*Li*Li/(4*self.par.wt1**2)

            self.max_temp = self.sim.max_wall_temp

            self.sigma = np.max(np.abs(np.array([temp_sigma_steady1,temp_sigma_steady2,temp_sigma_start])))
            if self.sigma == temp_sigma_start:
                self.max_temp = 293

            self.sigma_rat = self.met.sig_yield(self.max_temp)/(self.sigma*Kt(1+self.Di/self.par.wt1,self.par.rf2/self.par.wt1))
            return self.sigma_rat - self.sim.SF
        
        def pressure_drop_constraint(input):
            geom_update(input)
            flow_update()
            self.dp = deltaP(self.dhi,self.fi,self.fluid.rho,self.vi)
            return self.sim.pd_con - self.dp

        def mass_optimize(input):
            #print(input)
            geom_update(input)

            Lo =  np.linalg.norm(self.channelo.pointclasslist[0][0]-self.channelo.pointclasslist[1][0])
            def func(t):
                self.wto = t
                return np.sqrt(3/4*self.sim.p_cool*Lo*Lo/(1/self.sim.SF*self.met.sig_yield(self.sim.T_cool+100)*Kt(1+self.Do/self.wto,self.par.rf2/self.wto))) - self.wto
            scipy.optimize.root_scalar(func,bracket=[1e-4,1e-3])

            self.A_crosssection = np.pi*((self.par.ro+self.wto)**2 - (self.par.ri)**2) - (self.channeli.A + self.channelo.A) * self.N
            #print(self.A_crosssection)
            if pressure_drop_constraint(input) > 0 and stress_constraint(input) > 0:
                self.record_cross.append(self.A_crosssection)
                self.record_geom.append(input)
            return self.A_crosssection

        """
        nlc1 = scipy.optimize.NonlinearConstraint(stress_constraint,self.sim.SF,np.inf)
        nlc2 = scipy.optimize.NonlinearConstraint(pressure_drop_constraint,0,self.sim.pd_con)
        """
        cons = ({'type': 'ineq', 'fun': stress_constraint},
                {'type': 'ineq', 'fun': pressure_drop_constraint})
        x0 = np.array([self.par.t,self.par.wt1,self.par.rf1i,self.par.rf1o,self.par.rf2])
        #bounds =((1.5e-3,4e-3),(0.4e-3,2e-3),(0.1e-3,2e-3),(0.1e-3,2e-3),(0.15e-3,0.6e-3))
        bounds =((1.5e-3,4e-3),(0.4e-3,2e-3),(0.4e-3,4e-3),(0.1e-3,4e-3),(0.15e-3,0.6e-3))
        
        if pressure_drop_constraint(x0) > 0 and stress_constraint(x0) > 0:
            self.record_cross.append(mass_optimize(x0))
        try:
            scipy.optimize.minimize(mass_optimize,x0=x0,bounds=bounds,constraints=cons,method="SLSQP")
        except:
            input = self.record_geom[self.record_cross.index(min(self.record_cross))]
            geom_update(input)
            flow_update()
            pass

def geomi(x):
    cgp2 = np.array([0,x.ri+x.wt1])
    cgp1 = np.dot(ge.rot(-np.pi/x.N),np.array([0,x.ro-x.wt1]))
    cgp3 = np.dot(ge.rot(-2*np.pi/x.N),np.array([0,x.ri+x.wt1]))
    cgl1 = ge.line(cgp2,cgp1)
    cgl2 = ge.line(cgp1,cgp3)
    l1 = ge.shift(cgl1,x.wt2/2)
    l2 = ge.circle(x.ri+x.wt1)
    l3 = ge.shift(cgl2,x.wt2/2)
    f1 = ge.fillet_l_l(l1,l3,-1,-1,x.rf1i)
    f2 = ge.fillet_l_c(l1,l2,-1,-1,cgp2,x.rf2)
    f3 = ge.fillet_l_c(l3,l2,-1,-1,cgp3,x.rf2)


    p1 = ge.intersect(l1,ge.perp(l1,f1.p))
    p2 = ge.intersect(l1,ge.perp(l1,f2.p))
    p3 = l2.cord(ge.ang(cgp2,f2.p))
    p4 = l2.cord(ge.ang(cgp2,f3.p))
    p5 = ge.intersect(l3,ge.perp(l3,f3.p))
    p6 = ge.intersect(l3,ge.perp(l3,f1.p))

    shapelist= [[p1,l1],
                [p2,f2,1],
                [p3,l2,-1],
                [p4,f3,1],
                [p5,l3],
                [p6,f1,1]]
    
    return ge.shape(shapelist) , np.sqrt(np.linalg.norm(cgp1-p1)**2+x.rf1o**2)-x.rf1o

def geomo(x):
    cgp3 = np.array([0,x.ri+x.wt1])
    cgp1 = np.dot(ge.rot(-np.pi/x.N),np.array([0,x.ro-x.wt1]))
    cgp2 = np.dot(ge.rot(np.pi/x.N),np.array([0,x.ro-x.wt1]))
    cgl1 = ge.line(cgp2,cgp3)
    cgl2 = ge.line(cgp3,cgp1)
    l1 = ge.circle(x.ro-x.wt1)
    l2 = ge.shift(cgl1,-x.wt2/2)
    l3 = ge.shift(cgl2,-x.wt2/2)
    f1 = ge.fillet_l_c(l2,l1,1,1,cgp2,x.rf2)
    f2 = ge.fillet_l_l(l2,l3,1,1,x.rf1o)
    f3 = ge.fillet_l_c(l3,l1,1,1,cgp1,x.rf2)
    
    p1 = l1.cord(ge.ang(cgp3,f3.p))
    p2 = l1.cord(ge.ang(cgp3,f1.p))
    p3 = ge.intersect(l2,ge.perp(l2,f1.p))
    p4 = ge.intersect(l2,ge.perp(l2,f2.p))
    p5 = ge.intersect(l3,ge.perp(l3,f2.p))
    p6 = ge.intersect(l3,ge.perp(l3,f3.p))
    shapelist= [[p1,l1,1],
                [p2,f1,1],
                [p3,l2,],
                [p4,f2,1],
                [p5,l3],
                [p6,f3,1]]
    
    return ge.shape(shapelist), np.sqrt(np.linalg.norm(cgp1-p1)**2+x.rf1i**2)-x.rf1i


if __name__ == "__main__":

    import thermo

    ccase = sim(3000,400,5000,1.3e-6,70e5,50e5,0)
    in718 = metal(E=200e9,k=20,v=0.33,alpha=12e-6,sig_yield=([1150e6,1150e6,950e6,650e6,0],[273,50+273,700+273,850+273,1260+273]))

    t = 1.2e-3
    wt1 = 0.6e-3
    wt2 = 0.4e-3
    rf1 = 0.1e-3
    rf2 = 0.1e-3
    test = parameters(ri=65e-3,t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=48)

    fluid = thermo.Chemical('C2H5OH', T=288, P=50e5)

    t1  = time.time()
    test_res = physics(test,in718,ccase, fluid)
    t2 = time.time()
    print(t2-t1)

    print("f inner: " + str(test_res.fi))
    print("pressure drop: " + str(test_res.dp))
    print("velocity inner: " + str(test_res.vi))
    print("reynolds inner: " + str(test_res.Rei))
    print("hd inner: " + str(test_res.dhi))
    print("hd outer: " + str(test_res.dho))
    print("velocity outer: " + str(test_res.vo))
    print("reynold outer: " + str(test_res.Reo))
    print("ht inner: " + str(test_res.hi))
    print("ht outer: " + str(test_res.ho))
    print("T :" +str(test_res.par.t))
    print("wt1 : " +str(test_res.par.wt1))
    print("rf1i : " +str(test_res.par.rf1i))
    print("rf1o : " +str(test_res.par.rf1i))
    print("rf2 : " +str(test_res.par.rf2))
    print("wto : " +str(test_res.wto))