from pylab import *
from utils import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Solver(object):
    """
    
    Incompressible Navier Stokes Solver to solve Lid driven cavity.
    Written as the part of course CFD 1, AE 523

    Parameters:
    ----------------

    N = No. of internal grid points of pressure(!pressure)
    Re = Reynolds No.(default = 100.0)
    uw = Wall velocities <array of size 4> for all the four walls (default [0,1.0,0,0], upper wall moving)
    dt = time interval (default=1e-5)
    beta = Poisson solver parameter (default = 0.8)
    tol_p = Posson solver tolerance (default = 1e-4)
    tol_u = Velocity tolerance for final convergence check (default = 1e-6)
    maxiter = maximum number of iterations, in case convergence is not achieved(default = 50000)
    animate = if True, display dynamically chaging contour of u

    """
    def __init__(self, N):
        # N is number of pressure points in the staggered grid
        
        self.Ni = N + 2
        self.Nj = N + 2
        
        self.Li = 1.0
        self.Lj = 1.0
        self.dx = self.Lj/(N-1)
        self.dy = self.Li/(N-1)

        self.u = zeros([self.Ni+1, self.Nj], order="F")
        self.v = zeros([self.Ni, self.Nj-1], order="F")
        self.p = zeros([self.Ni, self.Nj], order="F")

        self.x = zeros([N, N], order="F")
        self.y = zeros([N, N], order="F")
        for i in range(N):
            for j in range(N):
                self.x[i,j] = i*self.dx
                self.y[i,j] = j*self.dy

        self.Re = 100.0
        self.dt = 1e-5
        self.beta = 0.8
        self.tol_p = 1e-4
        self.tol_u = 1e-6
        self.maxiter = 50000
        self.ut = np.zeros_like(self.u)
        self.vt = np.zeros_like(self.v)
        self.uw = np.zeros(4, order="F")
        self.uw[1] = 1.0
        self.p_old = np.zeros_like(self.p)
        self.u_old = copy(self.u)
        self.uu = np.zeros_like(self.x)
        self.vv = np.zeros_like(self.x)
        self.w = np.zeros_like(self.x)

    def check_convergence(self):
        self.residue = abs(self.u_old-self.u).max()
        if self.residue < self.tol_u and self.counter > 2 or self.counter > self.maxiter:
                return True
        else:
            self.counter +=1
            self.u_old = copy(self.u)
            if self.counter%100 == 0:
                print self.residue, self.counter
            return False

    def write_data(self):
        file_name = __file__.replace(".py","") + ".npz"
        np.savez(file_name, u=self.u, v=self.v, p=self.p, x=self.x, y=self.y, w=self.w)
        savetxt("w.txt", self.w)
        savetxt("v.txt", self.vv)
        savetxt("u.txt", self.uu)
        savetxt("x.txt", self.x)
        savetxt("y.txt", self.y)
        

        
    def step(self):
        set_boundary(self.u, self.v, self.uw)
        set_uv_t(self.u, self.v, self.ut, self.vt, self.Re, self.dx, self.dt)
        poisson(self.p, self.ut, self.vt, self.dx, self.dt, self.beta, self.tol_p)
        update_uv(self.u, self.v, self.ut, self.vt, self.p, self.dt, self.dx)
        calc_vorticity(self.u, self.v, self.uu, self.vv, self.w, self.dx, self.dy)

    def loop(self):
        self.counter = 0
        while 1:
            self.step()
            if self.check_convergence():
                self.write_data()
                break
            else:
                if self.counter%100 == 0:
                    self.write_data()
        

    def aniloop(self):
        self.counter = 0
        while 1:
            self.step()
            if self.check_convergence():
                self.write_data()
                break
            else:
                if self.counter%100 == 0:
                    self.write_data()
                    yield self.u[1:-1,1:-1],

        
    def run(self):
        if self.animate:
            fig, ax = plt.subplots()
            c = ax.quiver(self.x, self.y, self.uu, self.vv)
            def update(data):
                ax.cla()
                im = ax.quiver(self.x, self.y, self.uu, self.vv)
                title("Iteration: %i, Residue: %e"%(self.counter,self.residue))
                return im,
            
            ani = animation.FuncAnimation(fig, update, self.aniloop, interval=1)
            plt.show()
        else:
            self.loop()
            
if __name__ == "__main__":
    # init solver
    s = Solver(160)
    s.Re = 10.0
    s.uw = array([0.0,1.0,0.0,0.0])
    s.dt = .5e-5
    s.beta = 1.4
    s.tol_p = 1e-2
    s.tol_u = 1e-8
    s.maxiter = 500000
    s.animate = True
    s.run()

