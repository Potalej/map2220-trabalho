# def newton (p0, f, f_):
#     return p0 - f(p0)/f_(p0)

def newton (p0, f, f_p0):
    return p0 - f(p0)/f_p0

def derivada (f, x0):
    h = 0.001
    num = f(x0+h) - f(x0)
    den = h
    return num/den

from math import cos, sin, pi

f = lambda x: cos(x) - x
f_ = lambda x: -sin(x) - 1

p0 = pi/4

print('n \t pn')
print(0, '\t', p0)
for i in range(5):
    der = derivada(f, p0)
    p0 = newton(p0, f, der)
    print(i, '\t', p0)