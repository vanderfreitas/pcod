import matplotlib
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import matplotlib.figure
import math
import random
#from igraph import *


import sys

from matplotlib import rc
# activate latex text rendering
rc('text', usetex=True)



name = 'out_trajectory.csv'
csv_file = np.genfromtxt (name, delimiter="\t")




# parameters
N = 12
step = 1
beg = 0


iterations = int( (len(csv_file)/ float(N)) / step) - 1
print "iterations: ", iterations


#limites da tela (eixos x e y)
x_lim_inf = -45
x_lim_sup = 50
y_lim_inf = -20
y_lim_sup = 50


#x_lim_inf = -100
#x_lim_sup = 100
#y_lim_inf = -100
#y_lim_sup = 100


particulas = []
quivers = []



w, h = matplotlib.figure.figaspect(1.)
fig = pyplot.figure(figsize=(w,h))
ax = fig.add_subplot(111)
ax.set_xlim([x_lim_inf,x_lim_sup])
ax.set_ylim([y_lim_inf,y_lim_sup])

pyplot.xticks(fontsize=18)
pyplot.yticks(fontsize=18)


#ax.set_xlabel(r'$x$', usetex=True, fontsize=16)
#ax.set_ylabel(r'$y$', usetex=True, fontsize=16)



#Criando particulas
for i in xrange(N):
	index = beg*N*step + i 

	circle=pyplot.Circle((csv_file[index, 2],csv_file[index, 3]),2.0,fc='r', zorder=2)
	quivers.append(ax.quiver(csv_file[index, 2],csv_file[index, 3],0.001,0.001,angles='xy', scale_units='xy',scale=0.15, zorder=2))
	particulas.append(circle)




#Adicionando particulas na animacao
for i in xrange(N):
	index = beg*N*step + i
	particulas[i].center=(csv_file[index, 2],csv_file[index, 3])
	ax.add_patch(particulas[i])




R_x = 0
R_y = 0
old_R_x = 0
old_R_y = 0

for i in xrange(N):
    index = beg*N*step + i
    R_x = R_x + csv_file[index, 2]  
    R_y = R_y + csv_file[index, 3] 

R_x = R_x / float(N)
R_y = R_y / float(N)

old_R_x = R_x
old_R_y = R_y


#centro de massa
circle=pyplot.Circle((0,0),1.5,fc='b',zorder=2)
particulas.append(circle)
particulas[N].center=(R_x,R_y)
#ax.add_patch(particulas[N])



oldCoord = []
for i in xrange(N):
	index = beg*N*step + i
	oldCoord.append([csv_file[index, 2],csv_file[index, 3]])


def animate(i):
	global x_lim_inf
	global x_lim_sup
	global y_lim_inf
	global y_lim_sup
	global R_x
	global R_y
	global old_R_x
	global old_R_y
	global N

	print i



	old_R_x = R_x
	old_R_y = R_y

	R_x = 0.0
	R_y = 0.0



	for k in xrange(N):
		#####  Dynamics  #####
		index = (i+beg)*N*step + k

		old_x = oldCoord[k][0] #csv_file[index, 2]  #r_x[index];
		old_y = oldCoord[k][1] #csv_file[index, 3]  #r_y[index];


		R_x = R_x + old_x 
		R_y = R_y + old_y

		#Plotando o caminho percorrido pelas particulas, em cinza claro
		pyplot.plot([old_x, csv_file[index, 2]], [old_y, csv_file[index, 3]], '-', color=(0.7, 0.7, 0.7), zorder=0)


		particulas[k].center=(old_x,old_y)

		#Plotando os vetores de direcao das particulas
		quivers[k].set_offsets([old_x, old_y])
		quivers[k].set_UVC((csv_file[index, 2]-old_x),(csv_file[index, 3]-old_y))


		oldCoord[k][0] = csv_file[index, 2]
		oldCoord[k][1] = csv_file[index, 3]


	R_x = R_x / float(N)
	R_y = R_y / float(N)

	#pyplot.plot([old_R_x, R_x], [old_R_y, R_y], '-', color='blue', zorder=0)

	particulas[N].center=(R_x,R_y)
	#quivers[N].set_offsets([old_R_x, old_R_y])
	#quivers[N].set_UVC((R_x-old_R_x),(R_y-old_R_y))

	pyplot.tight_layout(pad=0.1)
	#pyplot.tight_layout(pad=0, w_pad=-5)
	#name = "fig_" + str(i) + ".pdf"
	name = '_fig_.png' # + str(i) + '.png'
	fig.savefig(name, dpi=100)


	return particulas



for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(25) 

for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(25)


ax.set_aspect(1)
pyplot.tight_layout(pad=0, w_pad=-5)



save_animation = False


if(save_animation):
	# Generating animation in mp4
	anim=animation.FuncAnimation(fig,animate,frames=iterations,blit=False)

	mywriter = animation.FFMpegWriter()
	anim.save('video.mp4',writer=mywriter)

else:
	# Iterating the data and saving it continuously as images
	for i in xrange(iterations):
		animate(i)




