
#import numpy as np
import mwlmc
Model = mwlmc.MWLMC()
X = Model.orbit([-8.27,0.,0.],[0.,240.,0.],1000,0.0005,127,127,127,True)
Model.print_orbit(X,'out.dat')

# test the different components
t,x,y,z = 0.,-8.27,0.,0.021
fx,fy,fz,dens,pot = Model.mwd_fields(t,x,y,z,127,False)
print('disc:',fx,fy,fz,dens,pot)
fx,fy,fz,dens,pot = Model.mwhalo_fields(t,x,y,z,127,False)
print('halo',fx,fy,fz,dens,pot)
fx,fy,fz,dens,pot = Model.lmc_fields(t,x,y,z,127,False)
print('lmc:',fx,fy,fz,dens,pot)
