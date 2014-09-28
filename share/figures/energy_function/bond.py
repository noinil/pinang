#!/usr/bin/env python

from pylab import *

figure(figsize=(8,5), dpi=80)
subplot(111)

X = np.linspace(0, 4, 256,endpoint=True)
C = 2*(X-2)*(X-2)

plot(X, C, color="green", linewidth=2.5, linestyle="-",
     label=r"$V_{bond} = K_r(r-r_0)^2$")

ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

xlim(-0.5, X.max()*1.1)
xticks([0, 2],
       [r'$0$', r'$r_0$'])

ylim(-0.8,5)
yticks([],
       [])

t = 2
plot([t,t],[0,4],
     color ='green',  linewidth=1., linestyle="--")

legend(loc='upper right')

savefig("bond.svg", dpi=72)
show()
