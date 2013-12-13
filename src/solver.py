from pylab import *
from utils import *
class Solver(object):
    """
    
    Incompressible Navier Stokes solver, mainly written to solve lid
    driven cavity problem, but can be modified to solve other problems as well.
    
    """
    def __init__(self):
        # N is number of pressure points in the staggered grid
        self.Ni = 80
        self.Nj = 80
        
        self.Li = 1.0
        self.Lj = 1.0

        self.dx = self.Lj/(self.Nj-1)
        self.dy = self.Li/(self.Ni-1)

        self.bctype = [0, 0, 0, 0]
        # 0: bottom face
        # 1: right face
        # 2: top face
        # 3: left face
        self.u = zeros([self.Ni+1, self.Nj], order="F")
        self.v = zeros([self.Ni, self.Nj-1], order="F")
        x = linspace(0, 1.0, self.Ni)
        y = linspace(0, 1.0, self.Ni)
        self.xx, self.yy = meshgrid(x, y)
        
        self.p = zeros([self.Ni, self.Nj], order="F")
        self.Re = 1000.0
        self.dt = 1e-5
        self.beta = 0.8
        self.tol = 1e-4
        self.ut = np.zeros_like(self.u)
        self.vt = np.zeros_like(self.v)
        self.uw = np.zeros(4, order="F")
        self.uw[1] = 1.0
        self.p_old = np.zeros_like(self.p)
        self.residue = 1.0

    def loop(self):
        u_old = copy(self.u)
        counter = 0
        while 1:
            set_boundary(self.u, self.v, self.uw)
            set_uv_t(self.u, self.v, self.ut, self.vt, self.Re, self.dx, self.dt)
            poisson(self.p, self.ut, self.vt, self.dx, self.dt, self.beta, self.tol)
            update_uv(self.u, self.v, self.ut, self.vt, self.p, self.dt, self.dx)
            self.residue = abs(u_old-self.u).max()
            if self.residue < 5e-6 and counter > 2 or counter > 10000:
                break
            else:
                u_old = copy(self.u)
                print self.residue, counter
                counter +=1
                if counter%1000 == 0:
                    np.savetxt("u.txt", self.u)
                    np.savetxt("v.txt", self.v)
                    np.savetxt("p.txt", self.p)
        contourf(self.u[1:-2,1:-2])
        colorbar()
        show()

if __name__ == "__main__":
    s = Solver()
    s.loop()
