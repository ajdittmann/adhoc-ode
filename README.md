### ADaptive-timestep High-Order Code for solving Ordinary Differential Equations 
# ADHOC ODE solver

Scipy's solve_ivp lacks some potentially useful features, like being able to set a custom minimum timestep as a function of state variables. This package provides similar basic ODE-solving functionality, with a few extra tools. It is a simple ad hoc ODE solver, and not intended to be particularly performant. More options and methods are to come, eventually.

## Installation
First, clone the repository using ``git clone https://github.com/ajdittmann/adhoc-ode.git``.  
Then, install the package using ``pip install -e .`` or ``python setup.py install``. 
