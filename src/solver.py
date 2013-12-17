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
        self.dx = self.Lj/(self.Nj-1)
        self.dy = self.Li/(self.Ni-1)

        self.u = zeros([self.Ni+1, self.Nj], order="F")
        self.v = zeros([self.Ni, self.Nj-1], order="F")
        self.p = zeros([self.Ni, self.Nj], order="F")

        x = linspace(0, 1.0, self.Ni)
        y = linspace(0, 1.0, self.Ni)
        self.xx, self.yy = meshgrid(x, y)
        for i in range(self.Ni):
            for j in range(self.Nj):
                self.xx[i,j] = i*self.dx
                self.yy[i,j] = j*self.dy

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
        np.savetxt("u.txt", self.u)
        np.savetxt("v.txt", self.v)
        np.savetxt("p.txt", self.p)
        np.savetxt("x.txt", self.xx)
        np.savetxt("y.txt", self.yy)
        
    def loop(self):
        self.counter = 0
        while 1:
            set_boundary(self.u, self.v, self.uw)
            set_uv_t(self.u, self.v, self.ut, self.vt, self.Re, self.dx, self.dt)
            poisson(self.p, self.ut, self.vt, self.dx, self.dt, self.beta, self.tol_p)
            update_uv(self.u, self.v, self.ut, self.vt, self.p, self.dt, self.dx)
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
            levels = linspace(-1.0,1.0,50)
            c = ax.contourf(self.u[1:-1,1:-1], levels)
            plt.colorbar(c)
            xlabel("x_index")
            ylabel("y_index")
            title("u-contour")
            def update(data):
                ax.cla()
                im = ax.contourf(np.rot90(data[0],k=3), levels)
                title("Iteration: %i, Residue: %e"%(self.counter,self.residue))
                return im,
            
            ani = animation.FuncAnimation(fig, update, self.loop, interval=1)
            #ani.save('u.mpeg', writer='ffmpeg', fps=10);
            plt.show()
        else:
            self.loop()

if __name__ == "__main__":
    # init solver
    s = Solver(80)
    s.Re = 1000.0
    s.uw = array([0.0,1.0,0.0,0.0])
    s.dt = .5e-7
    s.beta = 1.4
    s.tol_p = 1e-2
    s.tol_u = 1e-8
    s.maxiter = 500000
    s.animate = True
    s.run()

