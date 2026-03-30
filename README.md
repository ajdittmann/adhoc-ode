### ADaptive-timestep High-Order Code for solving Ordinary Differential Equations 
# ADHOC-ODE 

Scipy's solve_ivp lacks some potentially useful features, like being able to set a custom maximum stepsize as a function of state variables. This package provides similar basic ODE-solving functionality, with a few extra tools. This is a simple ad hoc ODE solver written in python, and not intended to be particularly performant. More options and methods are to come, eventually.

## Installation
First, clone the repository using ``git clone https://github.com/ajdittmann/adhoc-ode.git``.  
Then, install the package using ``pip install -e .``

## Usage
After installation, you can import the module using 
```
import adhocODE as adhoc

sol = adhoc.solve_ivp(fun, t_span, y0)
```
where ``fun`` is a function defined such that  ``dy/dt = f(t, y)``, ``t_span`` is a 2-element tuple with the starting and stopping times for the calculation, and ``y0`` is the initial value of the vector ``y``.

The full call signature of ``solve_ivp`` is
```
solve_ivp(fun, t_span, y0, args=None, tol=1e-8, t_eval=None, method="rk87", dtfunc=None, dtfunc_args=None)
```

Here, the optional arguments are
* ``args`` - *tuple*, additional arguments passed to ``fun``
* ``tol`` - *float*, target accuracy for ODE solution
* ``t_eval`` - *array-like* or *None*, times at which to output the approximate ODE solution
* ``method`` - *string*, which ODE solver to employ. Defaults to the high-order explicit Runge-Kutta method ``"rk87"`` [[1]](#references), but the 2nd-order implicit midpoint method ``"imid"`` is also supported
* ``dtfunc`` - *callable*, a function that returns a custom maximum timestep. This should have the same call signature as ``fun``, followed by any additional arguments
* ``dtfunc_args`` - *tuple*, additional arguments to ``dtfunc`` in addition to any passed using ``args``

## TO DO
* add other methods, particularly something implicit <!-- https://arxiv.org/pdf/2211.14574 looks promising  -->
* add option to check for going out-of bounds or NaNing, and retrying steps under certain conditions

## References
[[1]](https://epubs.siam.org/doi/10.1137/0715051) [[https://www.sfu.ca/~jverner/]](https://www.sfu.ca/~jverner/) J.H. Verner, SIAM NA 1978, 772-790,	"Explicit Runge-Kutta methods with estimates of the Local Truncation Error" 
