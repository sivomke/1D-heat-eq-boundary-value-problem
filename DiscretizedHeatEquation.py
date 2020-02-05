# defining heat process
# du_dt = 1/x^{m}*d_dx (x^{m}*du_dx) - q(x,t)*u + f(x,t), a<x<b, t>0
# m = 0: Cartesian coordinate system
# m = 1: Cylindrical coordinate system
# m = 2: Spherical coordinate system
# u = u(x, t) - temperature
# initial condition:
# u(x, 0) = u_0(x)
# boundary conditions:
# alpha_1*k(a, t)*du(a,t)_dx = beta_1*u(a, t) - mu_1(t)
# -alpha_2*k(b, t)*du(b,t)_dx = beta_2*u(b, t) - mu_2(t)

# sigma = 0 : explicit
# sigma = 0.5: Crank–Nicolson method
# sigma = 1: implicit

import numpy as np
import tridiagonal_matrix_algo as alg
import matplotlib.pyplot as plt

sigma = 0.5
R = 0.06  # cylinder radius
initial_temp = 293.15  # Kelvin
env_temp = 1093.15  # Kelvin
target_temp = 1073.15  # Kelvin ; x = 0 (at the cylinder center)
ro = 8900  # kg/m**3
thermal_cond = 380  # J/(kg*K)
gamma = 365  # W/(m**2*K)
K = 398  # W/(m*K) ; lambda
T_0 = 273.15  # 0 Celsius in Kelvin
N = 100  # 1/N - step, N+1 points
coord_points = np.linspace(0, 1, num=(N+1))
print(coord_points)
h = 1/N  # space step
# tau = 1/(N**2)  # time step
tau = 0.005
gamma_1 = gamma*R/K

alpha_1 = 1
beta_1 = 0
alpha_2 = 1
beta_2 = gamma_1
v_0 = 1

v_prev = np.array([v_0 for i in range(N+1)])
print(v_prev.size)


def c(i):
    if i == 0:
        return -alpha_1*0.5*h*(tau*sigma/h**2+0.5)
    elif i == N:
        return -sigma*tau*alpha_2*(1-0.5*h)/h**2-sigma*tau*beta_2/h-alpha_2*0.5*0.5*(2-h)
    else:
        return -(coord_points[i+1]+coord_points[i-1])*(0.5+tau*sigma/h**2)


def b(i):
    if i==0:
        return 0.5*tau*sigma*alpha_1/h
    else:
        return tau*sigma*(coord_points[i]+0.5*h)/h**2


def d(i):
    if i == N:
        return sigma*tau*alpha_2*(1-0.5*h)/h**2
    else:
        return tau*sigma*(coord_points[i]-0.5*h)/h**2


def phi(i, v_prev: np.ndarray):
    if i == 0:
        return -(1-sigma)*tau*alpha_1*0.5*h*(v_prev[1] - v_prev[0])/h**2 - 0.5*alpha_1*0.5*h*v_prev[0]
    elif i == N:
        return (1-sigma)*tau*alpha_2*(1-0.5*h)*(v_prev[N] - v_prev[N-1])/h**2 + tau*(1-sigma)*beta_2*v_prev[N]/h - \
                alpha_2*0.5*0.5*(2-h)*v_prev[N]
    else:
        return -0.5*(coord_points[i+1]+coord_points[i-1])*v_prev[i] -tau*sigma* \
               ((coord_points[i]+0.5*h)*(v_prev[i+1] - v_prev[i]) - (coord_points[i] - 0.5*h)*(v_prev[i]-v_prev[i-1]))/h**2


def p(i):
    return coord_points[i] - 0.5*h


trid_matrix = np.zeros((N+1,N+2))
trid_matrix[0,0] = c(0)
trid_matrix[0,1] = b(0)
trid_matrix[0, N+1] = phi(0, v_prev)
trid_matrix[N, N-1] = d(N)
trid_matrix[N,N] = c(N)
trid_matrix[N,N+1] = phi(N,v_prev)
for i in range(1, N):
    trid_matrix[i, i - 1] = d(i)
    trid_matrix[i, i] = c(i)
    trid_matrix[i, i+1] = b(i)
    trid_matrix[i, N + 1] = phi(i, v_prev)

u_prev = v_prev*(initial_temp-env_temp)+env_temp
print(u_prev)
iter = 0
while (u_prev[0]<=target_temp):
    trid_matrix[0, N + 1] = phi(0, v_prev)
    trid_matrix[N, N + 1] = phi(N, v_prev)
    for j in range(1, N):
        trid_matrix[j, N + 1] = phi(j, v_prev)
    v_prev = alg.tridAlgo(trid_matrix)
    u_prev = v_prev*(initial_temp-env_temp)+env_temp
    print(u_prev)
    iter +=1


print("iter: ",iter)
def plot():
    stop_time = tau*R*thermal_cond*ro/K
    x = coord_points*R
    func_val = v_prev*(initial_temp-env_temp)+env_temp - T_0
    plt.ylabel(u'Temp (\u00B0C)')
    plt.xlabel('Length  (metres)')
    plt.title("Temperature at {0:.2f} seconds ({1:.2f} hours) in the center = {2:.2f}".format(stop_time, stop_time/3600, func_val[0]))
    plt.plot(x, func_val)
    plt.show()


plot()












