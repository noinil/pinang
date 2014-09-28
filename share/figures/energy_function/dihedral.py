#!/usr/bin/env python

from pylab import *

figure(figsize=(8,5), dpi=80)
subplot(111)

X = np.linspace(-np.pi, 1.2*np.pi, 256,endpoint=True)
C = 1 - np.cos(X-np.pi/8)
S = 1 - np.cos( 3 * (X-np.pi/8))

plot(X, C, color="blue", linewidth=1.5, linestyle="--", label=r"$V_1 = 1-\cos(\phi-\phi_0)$")
plot(X, S, color="red", linewidth=1.5, linestyle="--",  label=r"$V_3 = 1-\cos 3(\phi-\phi_0)$")
plot(X, C+0.5*S, color="green", linewidth=2.5, linestyle="-",  label=r"$V_1 + 0.5 \times V_3$")

ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

xlim(X.min()*1.1, X.max()*1.1)
xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi, np.pi/8],
       [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$', r'$\phi_0$'])

ylim(-0.8,5)
yticks([+2, +1, +3],
       [r'$+2$', r'$+1$', r'$+3$'])

t = np.pi/8
plot([t,t],[0,3.2],
     color ='green',  linewidth=1., linestyle="--")

legend(loc='upper left')

savefig("dihedral.svg", dpi=72)
show()
