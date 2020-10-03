import numpy as np
import thermo
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

class LiquidInjector():
    def __init__(self, fluid, mixture, temperature, pressure, length, massflow, pressuredrop, inletangle):
        self.fluid = thermo.Mixture(fluid, ws=mixture,  T=temperature, P=pressure)
        self.length = length 
        self.massflow = massflow
        self.pressuredrop = pressuredrop
        self.inletangle = inletangle

    def xi1c(self, Re):
        # curfve fit from "Liquid Rocket Trhust Chambers" page 48 realting inlet efficiency to Reynolds number
        Re = np.log10(Re)
        return 3.55378*np.exp(-Re*0.647016) - 0.103358

    def xiinlet(self):
        return 0.5 + 1.2/np.pi*self.inletangle         # min 0.5 for coaxial flow before injection, max 0.9 for flow at pi/6 rad realtive to faceplate 
    
    def injector(self, maxiter=100, tol=1e-6):
        it = 0 
        difference = 1 
        xiinlet = self.xiinlet()
        xi = xiinlet
        while difference > tol:
            mu = 1/np.sqrt(1 + xi)
            diameter = 0.95*self.massflow**0.5 * mu**(-0.5) * (self.fluid.rho*self.pressuredrop)**(-0.25)
            vel = 1.273*self.massflow*self.fluid.rho**(-1)*diameter**(-2)
            Re = self.fluid.rho*vel*diameter/self.fluid.mu
            lam = 0.3164*Re**(-0.25)
            friction = lam*self.length/diameter
            xi1c = self.xi1c(Re)

            xi = xiinlet + xi1c + friction 
            newmu = 1/np.sqrt(1 + xi)

            it += 1
            if it > maxiter:
                raise ValueError("Not converged after ", maxiter, " iterations")
            difference = abs(mu - newmu)

        self.mu = mu
        self.diameter = diameter
        self.velocity = vel

class GasInjector():
    def __init__(self, gas, temperature, pressure, length, massflow, pressuredrop, upstreamdiameter, inletangle):
        self.gas = thermo.Chemical(gas, T=temperature, P=pressure)
        self.length = length 
        self.massflow = massflow
        self.pressuredrop = pressuredrop
        self.chamberpressure = pressure - pressuredrop
        self.upstreamdiameter = upstreamdiameter
        self.inletangle = inletangle

    def xiinlet(self):
        return 0.5 + 1.2/np.pi*self.inletangle
    
    # from "Liquid Rocket Trhust Chambers" page 52
    def injector(self, maxiter=100, tol=1e-6):
        if self.gas.phase == 'l':
            raise ValueError("Fluid is not gaseous before injection")
        it = 0 
        difference = 1 
        mu = 0.9                        # initial guess
        gamma = self.gas.Cpg/self.gas.Cvg
        R = self.gas.Cpg-self.gas.Cvg

        while difference > tol:
            lam2 = np.sqrt((gamma+1)/(gamma-1) * (1 - (self.chamberpressure/self.gas.P)**((gamma-1)/gamma)))
            q = ((gamma+1)/2)**(1/(gamma-1)) * lam2*(1 - (gamma-1)/(gamma+1)*lam2*lam2)**(1/(gamma-1))
            c = np.sqrt(gamma*R*self.gas.T) / (gamma*np.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1))))
            An = self.massflow*c / (mu*self.gas.P*q)
            diameter = 1.128*np.sqrt(self.massflow*c / (mu*self.gas.P*q))
            vel = lam2 * np.sqrt(2*gamma/(gamma-1) * R*self.gas.T * (1 - (self.chamberpressure/self.gas.P)**((gamma-1)/gamma)))
            Re = self.gas.rho*vel*diameter/self.gas.mu
            lam = 0.3164*Re**(-0.25)
            xi = lam*self.length/diameter + self.xiinlet()*(1-diameter**2/self.upstreamdiameter**2)
            newmu = 1/np.sqrt(1 + xi)

            a = np.sqrt(gamma*R*self.gas.T)
            if vel > a:
                self.pressuredrop -= 0.05e5
                print("flow velocity exceeding speed of sound, reducing design pressure drop to ", self.pressuredrop/1e6, " MPa")

            it += 1
            if it > maxiter:
                raise ValueError("Not converged after ", maxiter, " iterations")
            difference = abs(mu - newmu)
            mu = newmu

        self.mu = mu
        self.diameter = diameter
        self.velocity = vel


class AnnularOrifice():
    def __init__(self, fluid, temperature, pressure, length, pressuredrop, inner_radius, outer_radius):
        self.fluid = thermo.Chemical(fluid, T=temperature, P=pressure)
        self.length = length
        self.pressuredrop = pressuredrop
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius

    def injector(self):
        '''
        from https://www.mathworks.com/help/physmod/hydro/ref/annularorifice.html#bqvdf1z-1
        '''
        m = np.pi*self.outer_radius*(self.outer_radius - self.inner_radius)**3 / (6*self.fluid.rho*self.fluid.nu*self.length) * self.pressuredrop
        A = np.pi*(self.outer_radius**2 - self.inner_radius**2)
        mu = m / (A * np.sqrt(2*self.fluid.rho*self.pressuredrop))

        self.mu = mu
        self.massflow = m


class AnnulusInjector():
    def __init__(self, fluid, temperature, pressure, length, annulusdiameter, massflow, pressuredrop, inletangle):
        """UNVARIFIED model of annular injectors. Uses turbulent boundary layer theory to determine annular gap (di) and numerical solution to Darcey Weisbach eqaution to determine discharge coefficent 

        :param fluid: fluid flowing through annulus
        :type fluid: string
        :param length: length of annulus flow passage (coaxial distance)
        :param annulusdiameter: mean diameter of annulus
        :param pressuredrop: design pressure drop over injector 
        """        
        self.fluid = thermo.Chemical(fluid, T=temperature, P=pressure)
        self.annulusdiameter = annulusdiameter 
        self.length = length
        self.massflow = massflow
        self.pressuredrop = pressuredrop
        self.inletangle = inletangle

    def friction(self, Re, di):
        # solving Darcey Weisbach equation
        surface_roughness = 0
        fd = fsolve(lambda f: 1/(np.sqrt(f)) + 2*np.log10(surface_roughness/(3.7*di) + 2.51/(Re*np.sqrt(f))), 0.0000001)
        return fd[0]

    def xiinlet(self):
        #return 0.5 + 1.2/np.pi*self.inletangle             # min 0.5 for coaxial flow before injection, max 0.9 for flow at pi/6 rad realtive to faceplate 
        return 1.56                                          # fit from cryo test data 

    def injector(self, maxiter=100, tol=1e-6):
        it = 0 
        difference = 1 
        mu = 1/np.sqrt(1+self.xiinlet())                               
        di = 1e-3
        while difference > tol:
            D = self.annulusdiameter + di/2

            vel = mu*np.sqrt(2*self.pressuredrop/self.fluid.rho)
            Re = self.fluid.rho*vel*di/self.fluid.mu

            deltamax = 0*0.37*0.4*di/(Re**0.2)                # max thickness of turbulent boundary layer at 0.4 * annulus gap thickness
            di = self.massflow/(self.fluid.rho*np.pi*D*vel) + 1/8*deltamax

            xifriction = self.friction(Re, di)*self.length/di
            xiinlet = self.xiinlet()
            newmu = 1/np.sqrt(1+xifriction+xiinlet)

            it += 1
            if it > maxiter:
                raise ValueError("Not converged after ", maxiter, " iterations")
            difference = abs(mu - newmu)
            mu = newmu  

        #print(np.pi*D*di*newmu*np.sqrt(self.fluid.rho*2*self.pressuredrop))

        self.mu = mu
        self.diameter = di
        self.velocity = vel


class Pintle():
    def __init__(self, oxidiser_injector, fuel_injector):
        self.oxidiser_injector = oxidiser_injector
        self.fuel_injector = fuel_injector
        self.oxidiser_injector.injector()                  # compute injector properties 
        self.fuel_injector.injector()                      # compute injector properties 

    def momentum_ratio(self, oxidiser_angle, fuel_angle, n_oxidiser_holes, n_fuel_holes):
        axial_momentum = (self.fuel_injector.massflow * self.fuel_injector.velocity * np.sin(fuel_angle) * n_fuel_holes
                       + self.oxidiser_injector.massflow * self.oxidiser_injector.velocity * np.sin(oxidiser_angle) * n_oxidiser_holes
        )
        radial_momentum = (self.fuel_injector.massflow * self.fuel_injector.velocity * np.cos(fuel_angle) * n_fuel_holes
                       + self.oxidiser_injector.massflow * self.oxidiser_injector.velocity * np.cos(oxidiser_angle) * n_oxidiser_holes
        )

        self.tmr = radial_momentum/axial_momentum
        self.efficiency = (-2*self.tmr**2 + 1.4*self.tmr + 92)/100


def annulus_verification(inner_radius, outer_radius, fluid, pressure_range, temperature, discharge_coefficient):

    volumeflow = []
    D = 2*((outer_radius - inner_radius)/2 + inner_radius)
    di = outer_radius - inner_radius
    for p in pressure_range:
        liquid = thermo.Chemical(fluid, T=temperature, P=p)
        vel = discharge_coefficient*np.sqrt(2*p/liquid.rho)
        Re = liquid.rho*vel*di/liquid.mu
        deltamax = 0.37*0.4*di/(Re**0.2) 

        q = 2*np.pi*D*vel*(di/2 - 1/8*deltamax)
        q *= 60/0.001                               # conversion to l/min

        volumeflow.append(q)

    def experimental_volumeflow(pressuredrop):
        deltap = pressuredrop/1e5
        q_exp = (deltap/0.0012)**(1/2.1171)         # cryo test data in l/min
        return q_exp

    plt.plot(experimental_volumeflow(pressure_range), pressure_range/1e5, color = 'red', label = 'experimental volumeflow')
    plt.plot(volumeflow, pressure_range/1e5, color = 'blue', label = 'calculated volumeflow')
    plt.loglog()
    plt.grid()
    plt.xlabel('volume flow [l/min]')
    plt.ylabel('pressure drop [bar]')
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':

    n_holes = 4

    liq_inj = LiquidInjector(['h2o2', 'h2o'], [0.9,0.1], 288, 26.5e5, 3.46e-3, 0.163/n_holes, 6.5e5, np.pi/6)
    liq_inj.injector()
    print(liq_inj.mu)
    print(liq_inj.diameter*1000)

    gas_inj = GasInjector('o2', 288, 26e5, 4e-3, 0.163/4, 6e5, 20e-3, np.pi/2)
    gas_inj.injector()
    #print(gas_inj.mu)

    an_inj = AnnulusInjector('h2o', 288, 26.5e5, 2e-3, 24e-3/2, 0.447, 6.5e5, np.pi/2)
    an_inj.injector()
    #print(an_inj.mu)
    #print(an_inj.diameter*1000)

    annulus_verification(24.2e-3/2, 25.32e-3/2, 'h2o', np.arange(1e5, 15e5, 0.1e5), 288, 0.61)

    '''
    pressuredrops = np.arange(3e5, 10e5, 0.1e5)
    tmr = []
    for p in pressuredrops:
        liq_inj = LiquidInjector('o2', 90, 26.5e5, 3e-3, 0.667/n_holes, p)
        liq_inj.injector()
        pintle = Pintle(liq_inj, an_inj)
        pintle.momentum_ratio(0, np.pi/2, n_holes, 1)
        tmr.append(pintle.tmr)

    plt.plot(tmr, pressuredrops/1e5)
    plt.show()
    '''
