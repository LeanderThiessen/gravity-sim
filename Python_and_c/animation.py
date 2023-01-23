import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import sin, cos

class ball():
    def __init__(self):
        self.data = []

length = 10 #dimension of the system (square with sl length)
N=500        #number of balls
n=10         #only use every nth value in the data
def init():
    ax.set_ylim(-length, length)
    ax.set_xlim(-length, length)
    
print("Reading in data...")

fig, ax = plt.subplots()

balls = []
points = []


for i in range(N):
    balls.append(ball())

    data_file = np.genfromtxt("data\data_{}.txt".format(i))
    balls[i].data = data_file[::n]

    point, = ax.plot([0], [0], '.',markersize=3,color='red')
    point.set_data(0,0)
    points.append(point)


L= len(balls[0].data)
gen_list = []

for i in range(L):
    gen_list.append([balls[j].data[i] for j in range(N)])

gen_list = np.array(gen_list)

print("Start Simulation...")

def run(data):
    
    for i in range(N):
        x,y = data[i]
        points[i].set_data(x,y)

    

ani = animation.FuncAnimation(fig, run, gen_list, init_func=init, interval=1)
plt.show()