from model.ODES import ODES
from parallel_map import parallel_map
import numpy as np

odes1 = ODES ()
odes1.add_equation ("x", "x + y")
odes1.add_equation ("y", "x * y")
odes1.define_initial_value ("x", 1.0)
odes1.define_initial_value ("y", 1.0)

odes2 = ODES ()
odes2.add_equation ("x", "x / y")
odes2.add_equation ("y", "y / x")
odes2.define_initial_value ("x", 1.0)
odes2.define_initial_value ("y", 1.0)

t = np.array (range (20))
integrator = lambda sys : sys.evaluate_on (t)
parallel_map (integrator, ([odes1, odes2]) * 10, 4, pool_of_threads=True)
