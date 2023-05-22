import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from deffs import Runge_Kutta4, Addams

L = 4
m = 2

def g(t):
    return 9.81 + 0.01 * np.cos(2 * np.pi * t)

def gt(t):
    return -0.01 * np.sin(2*np.pi*t) * 2 * np.pi 

def right_part(t, s):
    x, y, vx, vy, T = s
    result = np.zeros_like(s)
    result[0] = vx
    result[1] = vy
    result[2] = -x/m/L*T
    result[3] = -y/m/L*T + g(t)
    result[4] = 2* m/L * (vx * result[2] + vy*result[3]) + m*gt(t)*y/L + m*g(t)*vy/L
    
    return result

t0, tf = 0, 4
iterations = 1000
vx = 4
state0 = np.array([0, L, vx, 0, m/L*(vx*vx) + m*g(0)])
step = (tf - t0)/iterations

begin_cond = Runge_Kutta4(right_part, 0, state0, 3, 0.01)

sol = solve_ivp(right_part, (t0, tf), state0, t_eval=np.linspace(t0, tf, 1000))

sol1 = Runge_Kutta4(right_part, 0, state0, iterations, step)

sol2 = Addams(right_part, 0, begin_cond, iterations, step)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(sol.y[0], sol.y[1])
ax.plot([item[0] for item in sol1], [item[1] for item in sol1], color='green')
ax.plot([item[0] for item in sol2], [item[1] for item in sol2], color='blue')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend(['y(x)'])
plt.show()


##result[4] = 2* m/L * (vx * result[2] + vy*result[3]) + m*gt(t)*y/L + m*g(t)*vy/L


##  T = m/L (vx^2 + vy^2) + mg(t)*cos(a)
