#!/bin/python

from pylab import *

figure(figsize=(8,5), dpi=80)
subplot(111)

X = np.linspace(0.5, 2, 256,endpoint=False)
C = 5*pow((1/X),12) - 6*pow((1/X),10)
# C1 = 5*pow((0.5/X),12) - 6*pow((0.5/X),10)
A = 5*pow((1/X),12)
R =  - 6*pow((1/X),10)

plot(X, C, color="green", linewidth=2.5, linestyle="-",
     label=r"$V_{native} = \varepsilon [5(\frac{\sigma}{r})^{12} - 6(\frac{\sigma}{r})^{10} ]$")
plot(X, A, color="blue", linewidth=1.5, linestyle="--",
     label=r"$V_{repulsion} = 5 \varepsilon (\frac{\sigma}{r})^{12}$")
plot(X, R, color="red", linewidth=1.5, linestyle="--",
     label=r"$V_{attraction} = - 6 \varepsilon (\frac{\sigma}{r})^{10}$")
# plot(X, C1, color="black", linewidth=1, linestyle="-.",
#      label=r"$V_{native} = \varepsilon [5(\frac{\sigma}{r})^{12} - 6(\frac{\sigma}{r})^{10} ]$")

ax = gca()
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.xaxis.set_ticks_position('top')
ax.spines['top'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0.8))

xlim(0.6, X.max()*1.01)
xticks([0.8,1],
       [r'$0.8\sigma$', r'$\sigma$'])

ylim(-3,5)
yticks([-1],
       [r'$-\varepsilon$'])

t = 1
plot([t,t],[0,-1],
     color ='green',  linewidth=1., linestyle="--")

plot([0.8,t],[-1,-1],
     color ='green',  linewidth=1., linestyle="--")

legend(loc='upper right')

savefig("native.svg", dpi=72)
show()
