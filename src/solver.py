class IncompressibleSolver(object):
    """
    
    Incompressible Navier Stokes solver, mainly written to solver lid
    driven cavity problem, but can be modified to solve other problems as well.
    
    """
    def __init__(self):
        self.Nx = 100
        self.Ny = 100
        self.Nb = 2
        
        self.Lx = 1.0
        self.Ly = 1.0

        self.dx = self.Lx/self.Nx
        self.dy = self.Ly/self.Ny

        # 0: bottom face
        # 1: right face
        # 2: top face
        # 3: left face
        self.u = zeros([self.Nx+2*self.Nb, self.Ny+2*self.Nb, 4])
        self.v = zeros([self.Nx+2*self.Nb, self.Ny+2*self.Nb, 4])
        self.p = ones([self.Nx+2*self.Nb, self.Ny+2*self.Nb])
        
    def apply_boundary_conditions(self):
        # bottom boundary
        for i in range(self.Nx+2*self.Nb):
            for j in range(self.Nb):
                if bctype == "solid": 
                    self.v[i,j,0] = - self.v[i,self.Nb,0]
                    self.u[i,j,0] = - self.u[i,self.Nb,0]
        # right boundary
        for i in range(self.Nx+self.Nb, self.Nx+2*self.Nb+1):
            for j in range(self.Ny+2*self.Nb):
                if bctype == "solid":
                    self.v[i,j,1] = - self.v[self.Nx+self.Nb-1,j,1]
                    self.u[i,j,1] = - self.u[self.Nx+self.Nb-1,j,1]
        # top boundary
        for i in range(self.Nx+2*self.Nb):
            for j in range(self.Ny+self.Nb, self.Ny+2*self.Nb+1):
                if bctype == "solid":
                    self.v[i,j,0] = - self.v[i,self.Ny+self.Nb-1,0]
                    self.u[i,j,0] = - self.u[i,self.Ny+self.Nb-1,0]
        # left boundary
        for i in range(self.Nb):
            for j in range(self.Ny+2*self.Nb):
                if bctype == "solid":
                    self.v[i,j,1] = - self.v[self.Nb,j,1]
                    self.u[i,j,1] = - self.u[self.Nb,j,1]

