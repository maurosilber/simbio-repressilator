import matplotlib.pyplot as plt
import numpy as np
from simbio.simulator import Simulator

from ...repressilator import Repressilator

t = np.arange(1000)
sim = Simulator(Repressilator(k=1, beta=2))
df = sim.run(t)
df.filter(like="protein").plot()
plt.show()
