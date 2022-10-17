""" In this program, we use the module shallow_water_simulation to create simulations
of the values outputted by the function shallow_water_simulation. The code saves
an animation that can be played. It also outputs the maximum time step that we can
use before our algorithm becomes unstable and checks whether the algorithm is stable
with the dt that we used.

We simulate a Gaussian peak centered at (x,y) = (0.5,0.5) in part (b) and a plane wave
in part (c)."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib.animation as anim
from mpl_toolkits.mplot3d import Axes3D
from shallow_water_simulation import check_stability, plus_and_minus, uvh, eta_star, shallow_water_simulation

""" === Question 1 === """
""" === Part b === """
Lx = 1 # length of spatial domain, x-direction in meters
Ly = 2 # length of spatial domain, y-direction in meters
T = 40 # length of time domain, in seconds
eps = 0.1 # smoothing parameter for eta


dx = 0.05 # grid spacing, x-direction
dy = 0.05 # grid spacing, y-direction
dt = 0.05 # grid spacing, in time

def eta_0(y,x):
    """Intial eta values function."""
    A = 0.05
    μ = 0.5
    σ = 0.049
    return A*np.exp(-1*((x-μ)**2 + (y-μ)**2)/(σ**2))

def H(y,x):
    """Bathymetry function."""
    return 0.001

# Grid Dimensions
N = int(Lx/dx) + 2 # number of grid cells, x-direction
M = int(Ly/dy) + 2 # number of grid cells, y-direction 
K = int(T/dt) + 1 # number of grid points, in time

# === Creating Initial Value Arrays ===
x, y = np.arange(0, N-1)*dx, np.arange(0, M-1)*dy
eta_init = np.zeros((M,N))
for i in range(1,M):
    eta_init[i,1:] = eta_0(np.ones(len(x))*y[i-1], x)
u_0 = np.zeros((M,N))
v_0 = np.zeros((M,N))
init_cond = [u_0, v_0, eta_init]
 
# === Stability ===
stable, dtmax = check_stability(dx, dy, dt, H, [0,Lx], [0,Ly], eta_init)

print("The largest possible value for dt is", np.round(dtmax,4))

if stable is True:
    print("The algorithm is stable.")
else:
    print("The algorithm is unstable.")
    
# === Calculating the Values for the Simulation ===
u,v,eta = shallow_water_simulation(Lx, Ly, T, eps, dx, dy, dt, init_cond, H)
    
# === Animation ===
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()  
ax = fig.add_subplot(111, projection='3d')

x, y = np.arange(1, N)*dx, np.arange(1, M)*dy
X,Y = np.meshgrid(x,y)
line = ax.plot_surface(X, Y, eta[0,1:,1:])

ax.set_xlim([0.0, M*dy])
ax.set_ylim([0.0, N*dx])
ax.set_zlim([-0.005, 0.02])
plt.title("Eta vs. (x,y)")
ax.set_xlabel("x-axis (m)")
ax.set_ylabel("y-axis (m)")
ax.set_zlabel("Eta (m)")

def animations(time):
    ax.clear()
    line = ax.plot_surface(X,Y,eta[time,1:,1:])
    ax.view_init(elev=30., azim=45)
    ax.set_xlabel("x-axis (m)")
    ax.set_ylabel("y-axis (m)")
    ax.set_zlabel("Eta (m)")
    ax.set_zlim([-0.005, 0.02])
    return line,

ani = FuncAnimation(fig, animations, frames=np.arange(0,K,10), interval=20000, blit=False)
ani.save('animation_part_b.mp4', writer=writer)
plt.show()

""" === Part c === """
eta_init = np.zeros((M,N))
eta_init[1,:] = 0.05
u_0 = np.zeros((M,N))
v_0 = np.zeros((M,N))
init_cond = [u_0, v_0, eta_init]

# === Stability ===
stable, dtmax = check_stability(dx, dy, dt, H, [0,Lx], [0,Ly], eta_init)

print("The largest possible value for dt is", np.round(dtmax,4))

if stable is True:
    print("The algorithm is stable.")
else:
    print("The algorithm is unstable.")
    
# === Calculating the Values for the Simulation ===
u,v,eta = shallow_water_simulation(Lx, Ly, T, eps, dx, dy, dt, init_cond, H)
  
#%%  
# === Animation ===
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()  
ax = fig.add_subplot(111, projection='3d')

x, y = np.arange(1, N)*dx, np.arange(1, M)*dy
X,Y = np.meshgrid(x,y)
line = ax.plot_surface(X, Y, eta[0,1:,1:])

ax.set_xlim([0.0, M*dy])
ax.set_ylim([0.0, N*dx])
ax.set_zlim([-0.005, 0.0])
plt.title("Eta vs. (x,y)")
ax.set_xlabel("x-axis (m)")
ax.set_ylabel("y-axis (m)")
ax.set_zlabel("Eta (m)")

def animations(time):
    ax.clear()
    line = ax.plot_surface(X,Y,eta[time,1:,1:])
    ax.view_init(elev=30., azim=45)
    ax.set_xlabel("x-axis (m)")
    ax.set_ylabel("y-axis (m)")
    ax.set_zlabel("Eta (m)")
    ax.set_zlim([-0.005, 0.15])
    return line,

ani = FuncAnimation(fig, animations, frames=np.arange(0,K,10), interval=20000, blit=False)
ani.save('animation_part_c.mp4', writer=writer)
plt.show()
    