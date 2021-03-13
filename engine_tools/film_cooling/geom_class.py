import numpy as np
from scipy import interpolate


def Kt(Dd,rd):
    #interpolation from https://web.njit.edu/~sengupta/met%20301/Stress%20Concentration.pdf trandeline in excel from 
    interp_const = np.array([[1.01,1.02,1.03,1.05,1.07,1.1,1.2,1.3,2,3,6],
                             [0.9327,0.9657,0.9896,0.9995,1.0149,1.0093,0.9884,0.9957,0.9811,0.9354,0.9141],
                             [-0.169,-0.179,-0.188,-0.199,-0.203,-0.217,-0.24,-0.256,-0.281,-0.318,-0.346]])

    A = interpolate.interp1d(interp_const[0],interp_const[1],fill_value='extrapolate')
    b = interpolate.interp1d(interp_const[0],interp_const[2],fill_value='extrapolate')
    return A(Dd)*rd**b(Dd)

class line:
    def __init__(self,p2,p1=np.array([0,0])):
        self.p2 = p2
        self.p1 = p1
        self.delta = (p2-p1)/np.linalg.norm(p2-p1)


class circle:
    def __init__(self,r1,pc= np.array([0,0])):
        self.r = r1
        self.p = pc

    def cord(self,angle):
        return np.dot(rot(angle),np.array([0,self.r])) + self.p
    
    def tangent(self,angle):
        vec = line(self.cord(angle),self.p)
        return perp(vec,self.cord(angle))


def perp(func,point):
    slope = np.dot(rot(np.pi/2),func.delta)
    return line(p2=slope+point,p1=point)

def ang(a,b):
    return np.arccos(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))*np.sign(a[0]*b[1]-a[1]*b[0])

def rot(angle):
    return np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])

def dis(func,distance):
    return func.delta*distance+func.p1

def shift(func,distance):
    ortho = perp(func,func.p1)
    return perp(ortho,dis(ortho,distance))

def intersect(line1,line2):
    M = np.array([[line1.delta[0],-line2.delta[0]],[line1.delta[1],-line2.delta[1]]])
    t1 = np.dot(np.linalg.inv(M),line2.p1-line1.p1)[0]
    return line1.delta*t1+line1.p1

def fillet_l_l(line1,line2,line1dir=1,line2dir=1,rf=0.1):
    line1new = shift(line1,line1dir*rf)
    line2new = shift(line2,line2dir*rf)
    return circle(rf,intersect(line1new,line2new))

def fillet_l_c(line_func,circle_func,linedir,circledir,ref,rf):
    angle = -ang(np.array([0,circle_func.r]),ref-circle_func.p)
    delta_angle = 1
    while np.abs(delta_angle*circle_func.r) > 10e-10:
        circletan = circle_func.tangent(angle)
        fillet = fillet_l_l(line_func,circletan,line1dir=linedir,line2dir=circledir,rf=rf)
        delta_angle = ang(circletan.p1,intersect(perp(circletan,fillet.p), circletan))*0.5
        angle += delta_angle
    return fillet

class shape:
    def __init__(self,pointclasslist):
        self.circ = 0
        self.A = 0
        for i in range(len(pointclasslist)):
            element = pointclasslist[i]
            nelement = pointclasslist[(i+1)%len(pointclasslist)]
            if isinstance(element[1],(line)):
                self.circ += np.linalg.norm(nelement[0] - element[0])
                self.A -= (nelement[0][0]-element[0][0])*((nelement[0][1]-element[0][1])/2 + element[0][1])
            
            elif isinstance(element[1],(circle)):
                angle = np.abs(ang(element[0]-element[1].p,nelement[0]-element[1].p))
                self.circ += angle*element[1].r
                self.A -= (nelement[0][0]-element[0][0])*((nelement[0][1]-element[0][1])/2 + element[0][1])
                self.A += element[1].r**2 * (angle/2 - np.sin(angle/2)*np.cos(angle/2))*element[2]
        self.pointclasslist = pointclasslist


class parameters:
    def __init__(self,ri,t,wt1,wt2,rf1,rf2,N):
        self.ri = ri
        self.t = t
        self.wt1 = wt1
        self.wt2 = wt2
        self.rf1i = rf1
        self.rf1o = rf1
        self.rf2 = rf2
        self.N = N
        self.ro = self.ri + self.t + self.wt1

    def local_channel_geometry(self):
        self.channeli , self.Do = geomi(self)
        self.channelo , self.Di = geomo(self)
        self.dhi = 4*self.channeli.A/self.channeli.circ
        self.dho = 4*self.channelo.A/self.channelo.circ

    def cooling_geometry(self,y):
        self.dhi_arr = np.ndarray(len(y))
        self.dho_arr = np.ndarray(len(y))
        self.Ai_arr = np.ndarray(len(y))
        self.Ao_arr = np.ndarray(len(y))
        self.wt1_arr = np.ones(len(y))*self.wt1
    
        for i in range(len(y)):
            self.ri = y[i]
            self.ro = y[i] + self.t + self.wt1
            self.local_channel_geometry()
            self.dhi_arr[i] = self.dhi
            self.dho_arr[i] = self.dho
            self.Ai_arr[i] = self.channeli.A
            self.Ao_arr[i] = self.channelo.A


class metal:
    def __init__(self,E,k,v,alpha,sig_yield):
        self.E = interpolate.interp1d(E[1],E[0],fill_value='extrapolate')
        self.k = k
        self.v = v
        self.alpha = alpha
        self.sig_yield_table = sig_yield
        self.sig_yield = interpolate.interp1d(sig_yield[1],sig_yield[0],fill_value='extrapolate')


class sim:
    def __init__(self, wall_temp, q, p_cool, p_chamber, SF=1.5):
        self.wall_temp = wall_temp
        self.q = q
        self.p_cool = p_cool
        self.p_cha = p_chamber
        self.SF = SF


class structural:
    def __init__(self,x,metal,sim):
        self.par = x
        self.sim = sim
        self.met = metal

    def geom_update(self):
        self.N = self.par.N
        self.channeli , self.Do = geomi(self.par)
        self.channelo , self.Di = geomo(self.par)
        self.dhi = 4*self.channeli.A/self.channeli.circ
        self.dho = 4*self.channelo.A/self.channelo.circ

    def stress(self):
        Li =  np.linalg.norm(self.channeli.pointclasslist[2][0]-self.channeli.pointclasslist[3][0])
        E = self.met.E(self.sim.wall_temp)

        temp_sigma_steady1 = E * self.met.alpha * self.sim.q * self.par.wt1 / (2*(1-self.met.v)*self.met.k) - (self.sim.p_cool-self.sim.p_cha)*Li*Li*3/(4*self.par.wt1**2)  
        temp_sigma_steady2 = E * self.met.alpha * self.sim.q * self.par.wt1 / (2*(1-self.met.v)*self.met.k) - (self.sim.p_cool-self.sim.p_cha)*Li*Li/(4*self.par.wt1**2)
        temp_sigma_start = (self.sim.p_cool)*Li*Li/(4*self.par.wt1**2)

        self.max_temp = self.sim.wall_temp

        self.sigma = np.max(np.abs(np.array([temp_sigma_steady1,temp_sigma_steady2,temp_sigma_start])))
        if self.sigma == temp_sigma_start:
            self.max_temp = 293

        self.sigma_rat = self.met.sig_yield(self.max_temp)/(self.sigma*Kt(1+self.Di/self.par.wt1,self.par.rf2/self.par.wt1))


def geomi(x):
    cgp2 = np.array([0,x.ri+x.wt1])
    cgp1 = np.dot(rot(-np.pi/x.N),np.array([0,x.ro-x.wt1]))
    cgp3 = np.dot(rot(-2*np.pi/x.N),np.array([0,x.ri+x.wt1]))
    cgl1 = line(cgp2,cgp1)
    cgl2 = line(cgp1,cgp3)
    l1 = shift(cgl1,x.wt2/2)
    l2 = circle(x.ri+x.wt1)
    l3 = shift(cgl2,x.wt2/2)
    f1 = fillet_l_l(l1,l3,-1,-1,x.rf1i)
    f2 = fillet_l_c(l1,l2,-1,-1,cgp2,x.rf2)
    f3 = fillet_l_c(l3,l2,-1,-1,cgp3,x.rf2)


    p1 = intersect(l1,perp(l1,f1.p))
    p2 = intersect(l1,perp(l1,f2.p))
    p3 = l2.cord(ang(cgp2,f2.p))
    p4 = l2.cord(ang(cgp2,f3.p))
    p5 = intersect(l3,perp(l3,f3.p))
    p6 = intersect(l3,perp(l3,f1.p))

    shapelist= [[p1,l1],
                [p2,f2,1],
                [p3,l2,-1],
                [p4,f3,1],
                [p5,l3],
                [p6,f1,1]]
    
    return shape(shapelist) , np.sqrt(np.linalg.norm(cgp1-p1)**2+x.rf1o**2)-x.rf1o

def geomo(x):
    cgp3 = np.array([0,x.ri+x.wt1])
    cgp1 = np.dot(rot(-np.pi/x.N),np.array([0,x.ro-x.wt1]))
    cgp2 = np.dot(rot(np.pi/x.N),np.array([0,x.ro-x.wt1]))
    cgl1 = line(cgp2,cgp3)
    cgl2 = line(cgp3,cgp1)
    l1 = circle(x.ro-x.wt1)
    l2 = shift(cgl1,-x.wt2/2)
    l3 = shift(cgl2,-x.wt2/2)
    f1 = fillet_l_c(l2,l1,1,1,cgp2,x.rf2)
    f2 = fillet_l_l(l2,l3,1,1,x.rf1o)
    f3 = fillet_l_c(l3,l1,1,1,cgp1,x.rf2)
    
    p1 = l1.cord(ang(cgp3,f3.p))
    p2 = l1.cord(ang(cgp3,f1.p))
    p3 = intersect(l2,perp(l2,f1.p))
    p4 = intersect(l2,perp(l2,f2.p))
    p5 = intersect(l3,perp(l3,f2.p))
    p6 = intersect(l3,perp(l3,f3.p))
    shapelist= [[p1,l1,1],
                [p2,f1,1],
                [p3,l2,],
                [p4,f2,1],
                [p5,l3],
                [p6,f3,1]]
    
    return shape(shapelist), np.sqrt(np.linalg.norm(cgp1-p1)**2+x.rf1i**2)-x.rf1i


if __name__ == '__main__':
    t = 2.4e-3
    wt1 = 0.55e-3
    wt2 = 0.6e-3
    rf1 = 0.3e-3
    rf2 = 0.4e-3
    test = parameters(ri=60e-3,t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=42)
    test.local_channel_geometry()
    print(test.dhi)
    print(test.dho)
    print(test.channeli.A)
    print(test.channelo.A)

    geometry = np.genfromtxt('sparrow_50bar.txt', delimiter='', dtype=None, skip_header = 13) / 1000
    test.cooling_geometry(geometry[:,1][::-1])
    print(test.dhi_arr)

    E_in718 = [199947961502.171,198569010043.535,195811107126.264,193053204208.992,190295301291.721,186847922645.132,184090019727.861,180642641081.271,177884738164,174437359517.411,170989980870.822,166853126494.915,163405747848.326,158579417743.101,153753087637.876,146858330344.698,139274097322.202,129621437111.752,119968776901.302,109626640961.535,98595029292.4497]
    T1_in718 = [294.3,310.9,366.5,422,477.6,533.2,588.7,644.3,699.8,755.4,810.9,866.5,922,977.6,1033.2,1088.7,1144.3,1199.8,1255.4,1310.9,1366.5]
    sig_in718 = [1123.85e6,1075.58e6,1020.42e6,965.27e6,930.79e6,799.79e6,689.48e6]
    T2_in718 = [293,588.7,810.9,922,977.6,1033.2,1088.7]
    in718 = metal(E=(E_in718,T1_in718),k=24,v=0.29,alpha=12e-6,sig_yield=(sig_in718, T2_in718))

    test = parameters(ri=25e-3,t=t,wt1=wt1,wt2=wt2,rf1=rf1,rf2=rf2,N=42)
    heat = sim(wall_temp=923, q=18e6*1.1, p_cool=72e5, p_chamber=28e5)
    struct = structural(test, in718, heat)
    struct.geom_update()
    struct.stress()
    print(struct.sigma_rat)