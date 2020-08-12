import numpy as np

class ptm_fields_3d(object):
    
    def __init__(self):
        
        return
    
    def set_grid(self,x,y,z):
    
        self.x = x
        self.y = y
        self.z = z
        self.nx = np.size(x)
        self.ny = np.size(y)
        self.nz = np.size(z)
        
        return
    
    def set_magnetic(self,bx,by,bz):
        
        bxs = np.shape(bx)
        bys = np.shape(by)
        bzs = np.shape(bz)
        
        if ~np.all([bxs==(self.nx,self.ny,self.nz),bys==bxs,bzs==bxs]):
            raise ValueError('dimension mismatch in set_magnetic')
        
        self.bx = bx
        self.by = by
        self.bz = bz
    
        return
    
    def set_electric(self,ex,ey,ez):
        
        exs = np.shape(ex)
        eys = np.shape(ey)
        ezs = np.shape(ez)
        
        if ~np.all([exs==(self.nx,self.ny,self.nz),eys==exs,ezs==exs]):
            raise ValueError('dimension mismatch in set_magnetic')
        
        self.ex = ex
        self.ey = ey
        self.ez = ez
        
        return
        
    def write_file(self,filename):
    
        with open(filename,'w') as f:
            f.write('{:4} {:4} {:4}\n'.format(self.nx,self.ny,self.nz))
            f.write((self.nx*'{:12.5e} ').format(*self.x).strip()+'\n')
            f.write((self.ny*'{:12.5e} ').format(*self.y).strip()+'\n')
            f.write((self.nz*'{:12.5e} ').format(*self.z).strip()+'\n')
            for i in range(self.nx):
                for j in range(self.ny):
                    for k in range(self.nz):
                        f.write(('{:4} {:4} {:4}'+(6*' {:12.5e}')).format(i+1,j+1,k+1,
                          self.bx[i,j,k],self.by[i,j,k],self.bz[i,j,k],
                          self.ex[i,j,k],self.ey[i,j,k],self.ez[i,j,k])+'\n')
        return
        
def binary_to_xyz(directory,runid):

    dir = directory if directory[-1]=='/' else directory+'/'
 
    xgrid=np.fromfile(dir+'xgrid.bin')
    ygrid=np.fromfile(dir+'ygrid.bin')
    zgrid=np.fromfile(dir+'zgrid.bin')
    
    nx,ny,nz=xgrid.size,ygrid.size,zgrid.size
    
    bx=np.fromfile(dir+'bx3d_{:04}.bin'.format(runid)).reshape([nx,ny,nz])
    by=np.fromfile(dir+'by3d_{:04}.bin'.format(runid)).reshape([nx,ny,nz])
    bz=np.fromfile(dir+'bz3d_{:04}.bin'.format(runid)).reshape([nx,ny,nz])
    ex=np.fromfile(dir+'ex3d_{:04}.bin'.format(runid)).reshape([nx,ny,nz])
    ey=np.fromfile(dir+'ey3d_{:04}.bin'.format(runid)).reshape([nx,ny,nz])
    ez=np.fromfile(dir+'ez3d_{:04}.bin'.format(runid)).reshape([nx,ny,nz])
    
    pf=ptm_fields_3d()
    pf.set_grid(xgrid,ygrid,zgrid)
    pf.set_magnetic(bx,by,bz)
    pf.set_electric(ex,ey,ez)
    pf.write_file(dir+'ptm_fields_{:04}.dat'.format(runid))
    
    return
    