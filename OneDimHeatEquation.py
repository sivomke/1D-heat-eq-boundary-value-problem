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


class OneDimHeatEquation:
    def __init__(self, m: int, k, q, f, a: float, b: float, u_0, alpha_1: int, alpha_2: int, beta_1: int, beta_2: int,
                 mu_1, mu_2):
        self.m = m
        self.k = k
        self.q = q
        self.f = f
        self.a = a
        self.b = b
        self.u_0 =  u_0
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.beta_1 = beta_1
        self.beta_2 = beta_2
        self.mu_1 = mu_1
        self.mu_2 = mu_2



