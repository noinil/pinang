#!/bin/python

from pylab import *

figure(figsize=(8,5), dpi=80)
subplot(111)

kk = 10.0
X = np.linspace(0, 2, 256,endpoint=False)
C = -1 + (1-exp(-kk*(X-1)))**2

A = kk*kk*(X-1)**2-1
B = X-X

plot(X, C, color="green", linewidth=2.5, linestyle="-",
     label=r"$V_{Morse} = -\varepsilon + \varepsilon(1-e^{-a(r-r_0)})^2$")
plot(X, A, color="blue", linewidth=2.5, linestyle="--",
     label=r"$V_{Morse} = -\varepsilon + \varepsilon a^2(r-r_0)^2$")

plot(X, B, color='cyan', linewidth=1.5, linestyle="--")

ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',-1))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0.8))

xlim(0.5, X.max()*1.01)
xticks([0.8,1],
       [r'$0.8r_0$', r'$r_0$'])

ylim(-1.2,0.5)
yticks([0,-1],
       [0,r'$-\varepsilon$'])

t = 1

plot([0.8,t],[-1,-1],
     color ='green',  linewidth=1., linestyle="--")

legend(loc='upper right')

savefig("morse.svg", dpi=72)
show()
