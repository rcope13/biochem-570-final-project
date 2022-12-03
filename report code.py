from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def SIIR_ODE(t, y):

    S,I1,I2,R1,R2,R3 = y # unpack y

    # parameters
    t1 = 4
    t2 = 10
    b = 0.4

    # reaction rates
    b1 = 0.3/t1
    b2 = 1/t1-b1
    c1 = 1/t2
    c2 = 0.8/t2
    c3 = 1/t2-c2
    d1 = 0.003
    
    # differential equations
    dSdt = -b*S*(I1+I2) + d1*R2
    dI1dt = b*S*(I1+I2) - (b1+b2)*I1
    dI2dt = b2*I1 - c1*I2
    dR1dt = b1*I1 - (c2+c3)*R1
    dR2dt = c1*I2 + c2*R1 - d1*R2
    dR3dt = c3*R1
    
    dydt = [dSdt, dI1dt, dI2dt, dR1dt, dR2dt, dR3dt] # repack dydt

    return dydt


# the time interval of the simulation
tspan = [0,100]

# set initial conditions
N = 1
I0 = 0.0004
y0 = [N-I0, I0, 0, 0, 0, 0] # S, I1, I2, R1, R2, R3

sol = solve_ivp(SIIR_ODE, tspan, y0)


plt.plot(sol.t,sol.y.T)
plt.legend(['S', 'I1', 'I2', 'R1', 'R2', 'R3'])
plt.xlabel('Time')
plt.ylabel('Proportion of Population')
plt.show()