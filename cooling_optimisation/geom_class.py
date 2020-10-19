import numpy as np
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