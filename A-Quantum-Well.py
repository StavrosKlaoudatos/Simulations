from matplotlib import animation, rc, pyplot as plt
from IPython.display import HTML
import numpy as np



N=1000
hbar=1
dt=0.1
m=1
xMin=0
xMax=500
x=np.linspace(xMin, xMax, N)
p=np.arange(-N/2, N/2)*((2*np.pi*hbar)/(N*((xMax-xMin)/(N-1))))
k0=m/2
x0=xMax/2
sigma=(xMax-xMin)/30
psi=(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-.5*((x-x0)/sigma)**2)*np.exp(1j*k0*x)


Vhat=np.zeros_like(x)
Vhat[350:499]=1
Vprop=np.exp(-1j*dt*Vhat/(2*hbar))
That=p**2/(2*m)
Tprop=np.exp(-1j*dt*That/hbar)

def Splitter(psi, Vprop=Vprop, Tprop=Tprop):
  psi*=Vprop
  psi_p=np.fft.fft(psi)
  psi_p*=Tprop
  psi=np.fft.ifft(psi_p)
  psi*=Vprop
  return psi





fig, ax = plt.subplots()
plt.close()

ax.set_xlim(( 0, 500))
ax.set_ylim((-0.001, 0.001))
line, = ax.plot([], [], lw=2)
lineV, = ax.plot([], [], lw=2)

def init():
  line.set_data([], [])
  lineV.set_data([], [])
  return (line, lineV,)
def animate(i):
  global psi
  psi=Splitter(psi)
  line.set_data(x, np.conj(psi)*psi)
  lineV.set_data(x, Vhat)
  return(line, lineV,)
anim = animation.FuncAnimation(fig, animate, init_func=init,
frames=400, interval=10, blit=True)

#For Google Colab
rc('animation', html='jshtml')
anim


