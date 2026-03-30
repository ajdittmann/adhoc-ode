import numpy as np

_METHODS = ["rk87", "imid"]

class IntegrationResult:
  """
  Container for outputs from the ODE solver.

  Attributes:
    t (np.ndarray): solution output times
    y (float): solution output values
    status (int): solver status
    message (str): output message
    success (bool): successful integration?
  """

  __slots__ = ("t", "y", "status", "message", "success")
  def __init__(self, t, y, status, message, success):
    self.t = t
    self.y = y
    self.status = status
    self.message = message
    self.success = success

def solve_ivp(fun, t_span, y0, args=None, tol=1e-8, t_eval=None, method="rk87", dtfunc=None, dtfunc_args=None):
  """Solve an initial value problem for a system of ODEs.

  This function numerically integrates a system of ordinary differential
  equations given an initial value::

    dy / dt = f(t, y)
    y(t0) = y0

  Here t is a 1-D independent variable (time), y(t) is an
  N-D vector-valued function (state), and an N-D
  vector-valued function f(t, y) determines the differential equations.
  The goal is to find y(t) approximately satisfying the differential
  equations, given an initial value y(t0)=y0.

  Some of the solvers support integration in the complex domain, but note
  that for stiff ODE solvers, the right-hand side must be
  complex-differentiable (satisfy Cauchy-Riemann equations [11]_).
  To solve a problem in the complex domain, pass y0 with a complex data type.
  Another option always available is to rewrite your problem for real and
  imaginary parts separately.

  Parameters
  ----------
  fun : callable
    Right-hand side of the system: the time derivative of the state ``y``
    at time ``t``. The calling signature is ``fun(t, y)``, where ``t`` is a
    scalar and ``y`` is an ndarray with ``len(y) = len(y0)``. Additional
    arguments need to be passed if ``args`` is used (see documentation of
    ``args`` argument). ``fun`` must return an array of the same shape as
    ``y``. See `vectorized` for more information.
  t_span : 2-member sequence
    Interval of integration (t0, tf). The solver starts with t=t0 and
    integrates until it reaches t=tf. Both t0 and tf must be floats
    or values interpretable by the float conversion function.
  y0 : array_like, shape (n,)
    Initial state. For problems in the complex domain, pass `y0` with a
    complex data type (even if the initial value is purely real).
  t_eval : array_like or None, optional
    Times at which to store the computed solution, must be sorted and lie
    within `t_span`. If None (default), use points selected by the solver.
  args : tuple, optional
    Additional arguments to pass to the user-defined functions.  If given,
    the additional arguments are passed to all user-defined functions.
    So if, for example, `fun` has the signature ``fun(t, y, a, b, c)``,
    then `jac` (if given) and any event functions must have the same
    signature, and `args` must be a tuple of length 3.
  tol : float or array_like, optional
    Absolute error tolerance, used to determine the step size.
  method : string, optional
    Which ODE solver to use. Currently only 'rk87' is supported.
  dtfunc : callable, optional
    A function that sets additional upper limits on the timestep. Its
    call signature should match that of `fun`, followed by additional
    optional arguments which can be passed using `dtfunc_args`.
  dtfunc_args : tuple, optional
    additional arguments for dtfunc, beyond those required by `fun` 

  Returns
  -------
  Bunch object with the following fields defined:
  t : ndarray, shape (n_points,)
    Time points.
  y : ndarray, shape (n, n_points)
    Values of the solution at `t`.
  status : int
    Reason for algorithm termination:

      * -1: Integration step failed.
      *  0: The solver successfully reached the end of `tspan`.
      *  1: A termination event occurred.

  message : string
    Human-readable description of the termination reason.
  success : bool
    True if the solver reached the interval end or a termination event
    occurred (``status >= 0``).

  References
  ----------
  .. [1] J.H. Verner, "Explicit Runge--Kutta methods with estimates of the
       Local Truncation Error", SIAM NA 1978, 772-790. 
  """

  ## check that valid ODE method was requested
  if method.lower() not in _METHODS:
    raise ValueError('method should be one of ', _METHODS)
  else:
    method = method.lower()

  t0, tf = map(float, t_span)
  tdir = np.sign(tf - t0)

  ## simplify passing args
  if args is not None:
    # Wrap the user's fun in lambdas to hide the additional
    # parameters.  Pass in the original fun as a keyword
    # argument to keep it in the scope of the lambda.
    try:
      _ = [*(args)]
    except TypeError as exp:
      suggestion_tuple = (
        "Supplied 'args' cannot be unpacked. Please supply `args`"
        f" as a tuple (e.g. `args=({args},)`)"
      )
      raise TypeError(suggestion_tuple) from exp

    def fun(t, x, fun=fun):
      return fun(t, x, *args)

    if dtfunc is not None:
      dt_args = list(args)
      if dtargs is not None: dt_args += dtfunc_args
      def dtfunc(t, x, dtfunc=dtfunc):
        return dtfunc(t, x, *dt_args)

  if args is None and dtfunc_args is not None:
    def dtfunc(t, x, dtfunc=dtfunc):
      return dtfunc(t, x, *dtfunc_args)


  ## check that t_eval is valid, if provided
  if t_eval is not None:
    t_eval = np.asarray(t_eval)
    if t_eval.ndim != 1:
      raise ValueError("`t_eval` must be 1-dimensional.")

    if np.any(t_eval < min(t0, tf)) or np.any(t_eval > max(t0, tf)):
      raise ValueError("Values in `t_eval` are not within `t_span`.")

    d = np.diff(t_eval)
    if tf > t0 and np.any(d <= 0) or tf < t0 and np.any(d >= 0):
      raise ValueError("Values in `t_eval` are not properly sorted.")

  ## pick ODE method
  if method == "rk87":
    from . import rk8
    solver = rk8.Solver(fun, len(y0) )
  else:
    from . import imid
    solver = imid.Solver(fun, len(y0) )

  ## pick try to pick initial timestep
  ## might need to make this robust...
  if dtfunc is not None:
    dt0 = dtfunc(t0, y0)
    dt0 *= tdir
  elif t_eval is not None:
    dt0 = (t_eval[1]-t_eval[0])*0.1
  else:
    dt0 = (t_span[1]-t_span[0])*0.001

  ## iterate a couple times to try and find a timestep yielding the desired tolerance
  f0, ee = solver.update(t0, y0, dt0)
  dt = solver.getDt(dt0, ee, tol)
  for i in range(2):
    f0, ee = solver.update(t0, y0, dt)
    dt = solver.getDt(dt, ee, tol)

  ts = []
  ys = []
  t_eval_i = 0

  tnow = t0
  ynow = y0
  if t_eval is None:
    t_eval = np.array([t0, tf])
    t_eval_i = 1

  else:
    if t0==t_eval[0]:
      ts.append(t0)
      ys.append(y0)
      t_eval_i = 1

  tnext = t_eval[t_eval_i]
  n_eval = len(t_eval)

  status = None
  while status is None:
    if (tnow + dt - tnext)*tdir >= 0:
      dt = tnext-tnow

      ynow, ee = solver.update(tnow, ynow, dt)
      tnow += dt
      t_eval_i += 1

      if t_eval_i <= n_eval:
        ts.append(tnow)
        ys.append(ynow)

      if t_eval_i < n_eval: tnext = t_eval[t_eval_i]
      else: tnext = tf
    else:
      ynow, ee = solver.update(tnow, ynow, dt)
      tnow += dt

    dt = solver.getDt(dt, ee, tol)
    if dtfunc is not None:
      dt = np.min([tdir*dt, dtfunc(tnow, ynow)])
      dt *= tdir

    if (tnow*tdir >= tf*tdir): status=0
    if np.isnan(dt): break

  ts = np.array(ts)
  ys = np.vstack(ys).T

  return IntegrationResult(ts, ys, status, "integration succeeded", 1)

