

#import numpy as np
import mwlmc
import matplotlib.pyplot as plt

Model = mwlmc.MWLMC()

X = Model.rewind((-8.27,0.,0.),(0.,240.,0.))
plt.figure()
plt.plot(X[0],X[1],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')



X = Model.mworbit((-8.27,0.,0.),(0.,200.,0.))
plt.figure()
plt.plot(X[0],X[2],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')


X = Model.get_lmc_trajectory()
plt.figure()
plt.plot(X[:,0],X[:,1],color='black',lw=1.)
plt.plot(X[:,0],X[:,2],color='black',lw=1.)
plt.plot(X[:,0],X[:,3],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')




#Model.get_expansion_centres_virial(0.0)

#X = Model.get_trajectories()
#plt.figure()
#plt.plot(X[:,0],X[:,1],color='black',lw=1.)
#plt.tight_layout()
#plt.savefig('orbittest.png')



X = Model.mworbit([-8.27,0.,0.],[0.,240.,0.],)
#Model.print_orbit(X,'out.dat')

plt.plot(X[0],X[1],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')


# test the different components
t,x,y,z = 0.,-8.27,0.,0.021
fx,fy,fz,dens,pot = Model.mwd_fields(t,x,y,z)
print('disc:',fx,fy,fz,dens,pot)
fx,fy,fz,dens,pot = Model.mwhalo_fields(t,x,y,z)
print('halo',fx,fy,fz,dens,pot)
fx,fy,fz,dens,pot = Model.lmc_fields(t,x,y,z)
print('lmc:',fx,fy,fz,dens,pot)

Model.get_expansion_centres_virial(t)
Model.get_expansion_centres_physical(t)
