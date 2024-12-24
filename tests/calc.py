from sympy import *

init_printing()
#EKF

q0,q1,q2,q3 = symbols("q0 q1 q2 q3")
wX,wY,wZ = symbols("wX wY wZ")
accX,accY,accZ = symbols("accX accY accZ")
dt = symbols("dt")



x = Matrix([q0,q1,q2,q3]) 

fx = Matrix([dt * 0.5 * (-wX * q1 - wY * q2 - wZ * q3), dt * 0.5 * (wX * q0 + wZ * q2 - wY * q3), dt * 0.5 * (wY * q0 - wZ * q1 + wX * q3), dt * 0.5 * (wZ * q0 + wY * q1 - wX * q2)])

hx = Matrix([2*(q1*q3-q0*q2),2*(q2*q3+q0*q1),1-2*(q1**2 +q2**2)])
z = Matrix([accX,accY,accZ])

print("\n")
pretty_print(Eq(x,fx))
print("\n")

print("\n")
pretty_print(fx.jacobian(x))
print("\n")

print("\n")
pretty_print(Eq(z,hx))
print("\n")

print("\n")
pretty_print(hx.jacobian(x))
print("\n")

