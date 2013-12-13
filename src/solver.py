from pylab import *
class Solver(object):
    """
    
    Incompressible Navier Stokes solver, mainly written to solve lid
    driven cavity problem, but can be modified to solve other problems as well.
    
    """
    def __init__(self):
        # N is number of pressure points in the staggered grid
        self.Ni = 40
        self.Nj = 40
        
        self.Nb = 2
        
        self.Li = 1.0
        self.Lj = 1.0

        self.dx = self.Lj/(self.Nj-1)
        self.dy = self.Li/(self.Ni-1)

        self.bctype = [0, 0, 0, 0]
        # 0: bottom face
        # 1: right face
        # 2: top face
        # 3: left face
        self.u = zeros([self.Ni+1, self.Nj])
        self.v = zeros([self.Ni, self.Nj-1])
        x = linspace(0, 1.0, self.Ni)
        y = linspace(0, 1.0, self.Ni)
        self.xx, self.yy = meshgrid(x, y)
        
        self.p = zeros([self.Ni, self.Nj])
        self.Re = 100.0
        self.dt = 1e-3
        self.beta = 0.8
    def apply_boundary_conditions(self):
        Ni = self.Ni
        Nj = self.Nj

        for i in range(Ni+1):
            for j in range(Nj):
                if j in [0]:
                    self.u[i,j] = -self.u[i,1]
                elif j in [Nj-1]:
                    self.u[i,j] = 1.0-self.u[i,Nj-2]
                if i in [0]:
                    self.u[i,j] = -self.u[1,j]
                elif i in [Ni+1-1]:
                    self.u[i,j] = -self.u[Ni+1-2,j]

        for i in range(Ni):
            for j in range(Nj-1):
                if j in [0]:
                    self.v[i,j] = -self.v[i,1]
                elif j in [Nj-1-1]:
                    self.v[i,j] = -self.v[i,Nj-1-2]
                if i in [0]:
                    self.v[i,j] = -self.v[1,j]
                elif i in [Ni-1]:
                    self.v[i,j] = -self.v[Ni-2,j] 
                    
    def calc_ut(self):
        self.ut = copy(self.u)
        for i in range(1, self.Ni+1-1):
            for j in range(1, self.Nj-1):
                a_1 = pow(0.5*(self.u[i+1,j] + self.u[i,j]), 2.0)
                a_2 = pow(0.5*(self.u[i,j] + self.u[i-1,j]), 2.0)
                a_3 = 0.25*(self.u[i,j] + self.u[i,j+1])*(self.v[i-1,j] + self.v[i,j])
                a_4 = 0.25*(self.u[i,j] + self.u[i,j-1])*(self.v[i-1,j-1] + self.v[i,j-1])

                d_1 = (self.u[i+1,j] - 2.0*self.u[i,j] + self.u[i-1,j])/self.dx/self.dx
                d_2 = (self.u[i,j+1] - 2.0*self.u[i,j] + self.u[i,j-1])/self.dy/self.dy

                A = (a_1 - a_2)/self.dx + (a_3 - a_4)/self.dy
                D = 1/self.Re*(d_1 + d_2)
                self.ut[i,j] = self.u[i,j] + self.dt*(-A+D)

        
    def calc_vt(self):
        self.vt = copy(self.v)
        for i in range(1, self.Ni-1):
            for j in range(1, self.Nj-1-1):
                a_1 = pow(0.5*(self.v[i,j+1] + self.v[i,j]),2.0)
                a_2 = pow(0.5*(self.v[i,j] + self.v[i,j-1]),2.0)
                a_3 = 0.25*(self.u[i+1,j] + self.u[i+1,j+1])*(self.v[i,j] + self.v[i+1,j])
                a_4 = 0.25*(self.u[i,j] + self.u[i,j+1])*(self.v[i-1,j] + self.v[i,j])

                d_1 = (self.v[i+1,j] - 2.0*self.v[i,j] + self.v[i-1,j])/self.dx/self.dx
                d_2 = (self.v[i,j+1] - 2.0*self.v[i,j] + self.v[i,j-1])/self.dy/self.dy
                
                A = (a_1 - a_2)/self.dy + (a_3 - a_4)/self.dx
                D = 1/self.Re*(d_1 + d_2)
                self.vt[i,j] = self.v[i,j] + self.dt*(-A+D)

    def calc_p(self):
        Ni = self.Ni
        Nj = self.Nj
        q = zeros_like(self.p)
        for i in range(1,Ni-1):
            for j in range(1,Nj-1):
                q[i,j] = ((self.ut[i+1,j] - self.ut[i,j])/self.dx + (self.vt[i,j] - self.vt[i,j-1])/self.dy)/self.dt
        p = copy(self.p)
        p_new = zeros_like(p)
        while 1:
            for i in range(1,Ni-1):
                for j in range(1,Nj-1):
                    p_new[i,j] = 0.25*self.beta*(p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1] - self.dx*self.dx*q[i,j]) + (1-self.beta)*p[i,j]

            if np.linalg.norm(p_new-p, np.inf) < 1e-3:
                for j in range(Nj):
                    p_new[0,j] = p_new[1,j]
                    p_new[Ni-1,j] = p_new[Ni-2,j]
                for i in range(Ni):
                    p_new[i,0] = p_new[i,1]
                    p_new[i,Nj-1] = p_new[i,Nj-2]

       
                p = copy(p_new)
                break
            else:
                p = copy(p_new)
        self.p = copy(p)

    def update_uv(self):
        for i in range(1, self.Ni+1-1):
            for j in range(1, self.Nj-1):
                self.u[i,j] = self.ut[i,j] - self.dt/self.dx*(self.p[i,j]-self.p[i-1,j])
        for i in range(1, self.Ni-1):
            for j in range(1, self.Nj-1-1):
                self.v[i,j] = self.vt[i,j] - self.dt/self.dy*(self.p[i,j+1]-self.p[i,j])
    
    def loop(self):
        u_old = copy(self.u)
        counter = 0
        while 1:
            self.apply_boundary_conditions()
            self.calc_ut()
            self.calc_vt()
            self.calc_p()
            self.update_uv()
            residue = abs(u_old-self.u).max()
            if residue < 1e-6:
                break
            else:
                print residue, counter
                u_old = copy(self.u)
                counter +=1
                if counter%10 == 0:
                    np.savetxt("u.txt", self.u)
                    np.savetxt("v.txt", self.v)
                    np.savetxt("p.txt", self.p)
        contourf(self.u[2:-3,2:-3])
        colorbar()
        show()

if __name__ == "__main__":
    s = Solver()
    s.loop()
