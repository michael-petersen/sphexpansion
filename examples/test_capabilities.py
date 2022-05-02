
# The harmonic flags are binary bit flags.
# That is, to turn on harmonic orders, add 2^(l-1) to the harmonic flag.
# For example, to enable l=1,2,3 to evolve, the harmonic flag should be set to 2^(1-1) + 2^(2-1) + 2^(3-1) = 7.
# To enable l=2,4 to evolve, the harmonic flag should be set to 2^(2-1) + 2^(4-1) = 10.
# To enable l=1,3,5, the harmonic flag should be set to 2^(1-1) + 2^(3-1) + 2^(5-1) = 21.
# The monopole is always on.

import mwlmc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

Model = mwlmc.MWLMC()

norbits = 10
pos = np.zeros([norbits,3])
vel = np.zeros([norbits,3])

for n in range(0,norbits):
    pos[n] = [1.*n,0.,0.]
    vel[n] = [0.,210.,0.]

X = Model.mworbit_parallel(pos,vel)

plt.figure()
for n in range(0,norbits):
    plt.plot(X[n][0],X[n][1],color=cm.viridis(n/(norbits-1.),1.),lw=1.)

plt.tight_layout()
plt.savefig('orbittest.png')



X = Model.rewind((-8.27,0.,0.),(0.,240.,0.),mwhharmonicflag=7,rewindtime=2.5,dt=0.001)
plt.figure()
plt.plot(X[0],X[1],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')


X = Model.mworbit((-8.27,0.,0.),(0.,200.,0.))
plt.figure()
plt.plot(X[0],X[1],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')


X = Model.get_lmc_trajectory()
plt.figure()
plt.plot(X[:,0],X[:,1],color='black',lw=1.)
plt.plot(X[:,0],X[:,2],color='black',lw=1.)
plt.plot(X[:,0],X[:,3],color='black',lw=1.)
plt.tight_layout()
plt.savefig('orbittest.png')
