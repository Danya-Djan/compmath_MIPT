import numpy as np

def Runge_Kutta2(Func, initialT, initialCondition, numOfIterations, s):
    solution = [initialCondition]
    t = initialT
    
    for i in range(numOfIterations):
        k1 = Func(t, solution[-1])
        k2 = Func(t + s, solution[-1] + s * k1)
        t += s
        solution.append(solution[-1] + s / 2 * (k1 + k2))
        
    return solution

def Runge_Kutta4(Func, initialT, initialСondition, numOfIterations, s):
    solution = [initialСondition]
    t = initialT
    
    for i in range(numOfIterations):
        k1 = Func(t, solution[-1])
        k2 = Func(t + s / 2, solution[-1] + s / 2 * k1)
        k3 = Func(t + s / 2, solution[-1] + s / 2 * k2)
        k4 = Func(t + s, solution[-1] + s * k3)
        t += s
        solution.append(solution[-1] + s / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
        
    return solution

def Addams(Func, initialT, initialСondition, numOfIterations, s):
    solution = initialСondition
    t = initialT
    
    for i in range(numOfIterations):
        solution.append(solution[-1] + s * (55/24*Func(t + 3*s, solution[-1]) - \
                        59/24*Func(t + 2*s, solution[-2]) + 37/24*Func(t + s, solution[-3]) - \
                        3/8*Func(t, solution[-4])))

        t += s
            
    return solution