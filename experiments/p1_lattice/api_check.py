import inspect, fpylll
from fpylll import BKZ
from fpylll.tools import bkz_simulator as s
print("fpylll", fpylll.__version__)
print("DEFAULT_STRATEGY", BKZ.DEFAULT_STRATEGY)
print([n for n in dir(s) if not n.startswith("_")])
print(inspect.signature(s.simulate))
print(inspect.getsource(s.simulate)[:2500])