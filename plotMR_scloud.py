import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Box1DKernel

data=np.genfromtxt('scloud0.txt', names=['P_z', 'T_z', 'Q_c', 'Q_v', 'Qt', 'Qs'])
P = data['P_z']
qc = data['Q_c']
qv = data['Q_v']
qs = data['Qs']
T = data['T_z']
qt = data['Qt']

data1=np.genfromtxt('scloud1.txt', names=['P_z', 'T_z', 'Q_c', 'Q_v', 'Qt', 'Qs'])
P1 = data1['P_z']
qc1 = data1['Q_c']
qv1 = data1['Q_v']
qs1 = data1['Qs']
T1 = data1['T_z']
qt1 = data1['Qt']

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
ax.semilogy(np.log10(qv),P,'k--', lw=3, label='Scloud = 0; Vapour')
ax.semilogy(np.log10(qc),P,'r+', lw=3, label='Scloud = 0; Condensate')
ax.semilogy(np.log10(qt),P,'m--', lw=3, label='Scloud = 0; Total')

#ax.semilogy(np.log10(qs1),P1,'c-', lw=3, label='Scloud = 1; Saturation')
ax.semilogy(np.log10(qv1),P1,'k-.', lw=3, label='Scloud = 1; Vapour')
ax.semilogy(np.log10(qc1),P1,'r.', lw=3, label='Scloud = 1; Condensate')
ax.semilogy(np.log10(qt1),P1,'b--', lw=3, label='Scloud = 1; Total')

ax.legend()
ax.set_xlim([-8, -4])
ax.set_ylim([1, .1])
#ax.set_xlabel('log10(Volume mixing ratio)', fontsize=fontsize)
#ax.set_ylabel('Pressure [bar]', fontsize=fontsize)

#ax2 = ax.twiny()

#ax.plot(np.log10(qc_smooth),P,'k--')

#ax2.set_xlim([110, 160])
fig.savefig('mixingR_Scloud.png',bbox_inches='tight')

plt.show()
