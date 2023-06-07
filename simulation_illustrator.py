import numpy as np
import os
import matplotlib.pyplot as plt

path = "C:\\Science\\1M"
filename = 'system_state704.npz'
os.chdir(path)
state = np.load(filename)

x = state['x'].transpose()
y = state['y'].transpose()

plt.figure()
plt.title('System states')
plt.xlabel('AU')
plt.ylabel('AU')
legend = f"Modelling time: {state['t'][-1]} years"
plt.plot(x[0,:], y[0,:], label = legend)
plt.legend(loc = 'best', fontsize = 7)
for i in range(1, len(state['m'])):
    plt.plot(x[i,:], y[i,:])
os.chdir('C:\\Users\\nmago\\Desktop')
plt.xlim(-20, 60)
plt.ylim(-40, 40)
plt.savefig('flyby.pdf')
plt.show()

plt.figure()
plt.title('Energy changing', fontsize = 16)
plt.xlabel('Time, years', fontsize = 14)
plt.ylabel(r'Energy,  $\frac{AU^2 * M_e}{d^2}$', fontsize = 14)
K = state['K'].transpose()
P = state['P'].transpose()
t = state['t']
for i in range(len(state['m']) - 1):
    plt.plot(t, K[i,:] + P[i,:])
K_sum, P_sum = 0, 0
for i in range(len(state['m']) - 1):
    K_sum += K[i]
    P_sum += P[i]
plt.plot(t, K_sum + P_sum)
plt.legend((f"first, m = {state['m'][0]}", f"second m = {state['m'][1]}", 'total'))
plt.show()
