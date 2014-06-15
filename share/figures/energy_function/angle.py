#!/bin/python

from pylab import *

figure(figsize=(8,5), dpi=80)
subplot(111)

X = np.linspace(-2, 4, 256,endpoint=True)
C = 0.5*(X-1)*(X-1)

plot(X, C, color="green", linewidth=2.5, linestyle="-",
     label=r"$V_{angle} = K_\theta(\theta-\theta_0)^2$")

ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

xlim(X.min()*1.2, X.max()*1.2)
xticks([-2, 0, 2, 4, 1],
       [r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$', r'$\theta_0$'])

ylim(-0.8,5)
yticks([],
       [])

t = 1
plot([t,t],[0,4],
     color ='green',  linewidth=1., linestyle="--")

legend(loc='upper right')

savefig("angle.svg", dpi=72)
show()
