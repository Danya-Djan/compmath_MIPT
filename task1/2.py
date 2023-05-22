import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from deffs import Runge_Kutta2, Runge_Kutta4

# Объявление констант
# m = 0.01215088
m = 0.012277471
M = 1 - m
f = 1

# Объвление вспомогательных функций
def r_1(x, y):
    return np.sqrt((x + m)**2 + y**2)

def r_2(x, y):
    return np.sqrt((x - M)**2 + y**2)

# Функция правой части дифференциального уравнения y' = f(t, y)
def rightPart(t, s):
    x, y, vx, vy = s
    # dx
    dx = vx
    # dy
    dy = vy
    # dvx
    dvx = 2*vy + x - M*(x+m)/r_1(x, y)**3 - m*(x - M)/r_2(x, y)**3 - f*vx
    # dvy
    dvy = -2*vx + y - M*y/r_1(x, y)**3 - m*y/r_2(x, y)**3 - f*vy
    return np.array([dx, dy, dvx, dvy])


# Проинициализируем значение шага и количество иттераций
iterations = 30000
beg_t = 0
end_t = 16
step = (end_t - beg_t)/iterations

# Проинициализируем начальные значения
x_0 = 1.2
y_0 = 0
dx_0 = 0
dy_0 = -1.05
# x_0 = 0.994
# y_0 = 0
# dx_0 = 0
# dy_0 = -2.031732629557337
state0 = np.array([x_0, y_0, dx_0, dy_0])

# Вызов метода численного решения
t_span = (beg_t, end_t)
sol = solve_ivp(rightPart, t_span, state0, method='RK45', rtol=1e-12, atol=1e-12, t_eval=np.linspace(beg_t, end_t, iterations))
sol1 = Runge_Kutta4(rightPart, 0, state0, iterations, step)
sol2 = Runge_Kutta2(rightPart, 0, state0, iterations, step)



# Plotting the result
fig, ax = plt.subplots(figsize=(8, 6))
plt.title('Trajectory of the restricted three-body problem')
ax.plot(sol.y[0], sol.y[1], color='red')
#ax.plot(sol1.y[0], sol1.y[1], color='green')
ax.plot([item[0] for item in sol1], [item[1] for item in sol1], color='blue')
ax.plot([item[0] for item in sol2], [item[1] for item in sol2], color='green')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.axis('equal')
plt.show()
