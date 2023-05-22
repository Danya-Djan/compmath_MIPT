import numpy as np
import matplotlib.pyplot as plt

time_periods = [1, 80, 160, 240, 320, 400]

def u_0x(x):
    if x >= 1 and x <= 3:
        return 1
    else:
        return 0

def u_t0(t):
    return 0

K = 400
N = 400
T = 8
X = 10
tau = T / K
h = X / N

t = np.linspace(0, T, K + 1)
x = np.linspace(0, X, N + 1)

u = np.zeros((K + 1, N + 1))
for n in range(N + 1):
    u[0, n] = u_0x(x[n])
for k in range(K + 1):
    u[k, 0] = u_t0(t[k])


# Схема Лакса
for k in range(K):
    for n in range(1, N):
        u[k + 1, n] = 0.5 * (u[k, n + 1] + u[k, n - 1]) / (0.5 * tau * (u[k, n + 1] - u[k, n - 1]) / h + 1)

plt.figure(figsize=(12, 8))
plt.rc('font', **{'size': 14})
plt.suptitle('Схема Лакса')



for i, k in enumerate(time_periods):
    plt.plot(x, u[k, :], label=f'Time Step {k}')

plt.xlim([0, 10])
plt.xticks(np.linspace(0, 10, 6))
plt.xlabel('$x$')
plt.ylim([0, 2])
plt.yticks(np.linspace(0, 2, 6))
plt.ylabel('$u$')
plt.grid()
plt.legend()
plt.savefig('pics/Схема Лакса.png')
plt.close()


# Схема Лакса-Вендроффа
for k in range(K):
    for n in range(1, N):
        z = u[k, n]
        for s in range(7):
            z = z - (z - u[k, n] + 0.5 * tau * z * \
                     (u[k, n + 1] - u[k, n - 1] - tau * z * (u[k, n + 1] - 2 * u[k, n] + u[k, n - 1]) / h) / h) / \
            (1 + 0.5 * tau * (u[k, n + 1] - u[k, n - 1]) / h - \
             tau**2 * z * (u[k, n + 1] - 2 * u[k, n] + u[k, n - 1]) / h**2)
        u[k + 1, n] = z

plt.figure(figsize=(12, 8))
plt.rc('font', **{'size': 14})
plt.suptitle('Схема Лакса-Вендроффа')



for i, k in enumerate(time_periods):
    plt.plot(x, u[k, :], label=f'Time Step {k}')

plt.xlim([0, 10])
plt.xticks(np.linspace(0, 10, 6))
plt.xlabel('$x$')
plt.ylim([0, 2])
plt.yticks(np.linspace(0, 2, 6))
plt.ylabel('$u$')
plt.grid()
plt.legend()
plt.savefig('pics/Схема Лакса-Вендроффа.png')
plt.close()


# Явный левый уголок
for k in range(K):
    for n in range(0, N -1):
        u[k + 1, n] = u[k, n] * (1 - tau * (u[k, n + 1] - u[k, n ]) / h)

plt.figure(figsize=(12, 8))
plt.rc('font', **{'size': 14})
plt.suptitle('Явный левый уголок')



for i, k in enumerate(time_periods):
    plt.plot(x, u[k, :], label=f'Time Step {k}')

plt.xlim([0, 10])
plt.xticks(np.linspace(0, 10, 6))
plt.xlabel('$x$')
plt.ylim([0, 2])
plt.yticks(np.linspace(0, 2, 6))
plt.ylabel('$u$')
plt.grid()
plt.legend()
plt.savefig('pics/Явный левый уголок.png')
plt.close()


# Явный правый уголок
for k in range(K):
    for n in range(1, N):
        u[k + 1, n] = u[k, n] / (1 + tau * (u[k, n + 1] - u[k, n]) / h)

plt.figure(figsize=(12, 8))
plt.rc('font', **{'size': 14})
plt.suptitle('Явный правый уголок')



for i, k in enumerate(time_periods):
    plt.plot(x, u[k, :], label=f'Time Step {k}')

plt.xlim([0, 10])
plt.xticks(np.linspace(0, 10, 6))
plt.xlabel('$x$')
plt.ylim([0, 2])
plt.yticks(np.linspace(0, 2, 6))
plt.ylabel('$u$')
plt.grid()
plt.legend()
plt.savefig('pics/Явный правый уголок.png')
plt.close()