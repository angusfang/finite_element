import numpy as np
import matplotlib.pyplot as plt

distance=np.linspace(0,20,500)
stress_maxs=np.load('stress_maxs.npy')
plt.plot(distance,stress_maxs)
print(np.min(stress_maxs))
plt.show()