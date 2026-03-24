### ADaptive-timestep High-Order Code for solving Ordinary Differential Equations 
# ADHOC ODE solver

Scipy's solve_ivp lacks some potentially useful features, like being able to set a custom minimum timestep as a function of state variables. This package provides similar basic ODE-solving functionality, with a few extra tools. It is a simple ad hoc ODE solver, and not intended to be particularly performant. More options and methods are to come, eventually.

## Installation
First, clone the repository using ``git clone https://github.com/ajdittmann/adhoc-ode.git``.  
Then, install the package using ``pip install -e .`` or ``python setup.py install``. 

## Usage
After installation, you can import the module using 
```
import adhocODE as adhoc

sol = adhoc.solve_ivp(fun, t_span, y0)
```
where ``fun`` is a function defined such that  ``dy/dt = f(t, y)``, ``t_span`` is a 2-element tuple with the stopping and starting times for the simulation, and ``y0`` is the initial value of the vector ``y``.

The full call signature of ``solve_ivp`` is
```
solve_ivp(fun, t_span, y0, args=None, tol=1e-8, t_eval=None, dtfunc=None, method="rk87")
```

Here, the optional arguments are
* ``args``
* ``tol``
* ``t_eval``
* ``dtfunc``
* ``method``

## TO DO
* add other methods, particularly something implicit
* add option to check for going out-of bounds or NaNing, and retrying steps under certain conditions
