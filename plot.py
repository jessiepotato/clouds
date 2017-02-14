import matplotlib.pyplot as plt
import numpy as np

data=np.genfromtxt('output.txt', names=['P_z', 'x'])
P = data['P_z']
x = data['x']


fig = plt.figure()

ax = fig.add_subplot(111)

ax.invert_yaxis()

#ax.plot(Qc,P,'.')

#ax.plot(T,P,'r.')

#qc_smooth = convolve(qc, Box1DKernel(10))
#qc_0 = qc_smooth==0
#qc_smooth[qc_0]=10E-10

#qc_0 = qc==0
#qc[qc_0]=10E-10

#qv_0 = qv==0
#qv[qv_0]=10E-10

fontsize=25
#ax.semilogy(np.log10(qs),P,'c-', lw=3, label='Scloud = 0; Saturation')
ax.semilogy(x,P,'k.-', lw=3, label='')


#ax.legend()
#ax.set_xlim([-8, -4])
#ax.set_ylim([1, .1])
ax.set_xlabel('Geometric mean radius [micron]', fontsize=fontsize)
#ax.set_xlabel('alpha', fontsize=fontsize)
ax.set_ylabel('Pressure [bar]', fontsize=fontsize)

#ax2 = ax.twiny()

#ax.plot(np.log10(qc_smooth),P,'k--')

#ax2.set_xlim([110, 160])
fig.savefig('mixingR_Scloud.png',bbox_inches='tight')

plt.show()
