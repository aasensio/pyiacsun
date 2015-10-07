__all__ = ['deriv']

import numpy as np

def deriv(xIn, y):
	"""
	The DERIV function uses three-point (quadratic) Lagrangian interpolation to compute the derivative of an evenly-spaced or unevenly-spaced array of data.
	Given three neighboring points in your data with X coordinates [x0, x1, x2] and Y coordinates [y0, y1, y2], the three-point Lagrangian interpolation polynomial is defined as:
	y = y0(x – x1)(x – x2)/[(x0 – x1)(x0 – x2)] +
   	y1(x – x0)(x – x2)/[(x1 – x0)(x1 – x2)] +
  		y2(x – x0)(x – x1)/[(x2 – x0)(x2 – x1)]
	This can be rewritten as:
		y = y0(x – x1)(x – x2)/(x01x02) – y1(x – x0)(x – x2)/(x01x12) + y2(x – x0)(x – x1)/(x02x12)
		where x01 = x0 – x1, x02 = x0 – x2, and x12 = x1 – x2.
	Taking the derivative with respect to x:
		y' = y0(2x – x1 – x2)/(x01x02) – y1(2x – x0 – x2)/(x01x12) + y2(2x – x0 – x1)/(x02x12)

	Given a discrete set of X locations and Y values, the DERIV function then computes the derivative at all of the X locations. For example, for all of the X locations (except the first and last points), the derivative y' is computed by substituting in x = x1:
		y' = y0x12/(x01x02) + y1(1/x12 – 1/x01) – y2x01/(x02x12)
	The first point is computed at x = x0:
		y' = y0(x01 + x02)/(x01x02) – y1x02/(x01x12) + y2x01/(x02x12)
	The last point is computed at x = x2:
		y' = –y0x12/(x01x02) + y1x02/(x01x12) – y2(x02 + x12)/(x02x12)
	
	Args:
	    xIn (TYPE): x-axis
	    y (TYPE): y-axis
	
	Returns:
	    TYPE: derivatives computed along the first axis
	"""
	shape = y.shape
	ndim = len(shape)-1

	x = np.copy(xIn)
	for i in range(ndim):
		x = np.expand_dims(x, axis=1)

	n = shape[0]
	d = np.zeros(shape)

	x12 = x - np.roll(x, -1, axis=0)
	x01 = np.roll(x, 1, axis=0) - x
	x02 = np.roll(x, 1, axis=0) - np.roll(x, -1, axis=0)

	d = np.roll(y, 1, axis=0) * (x12 / (x01*x02)) + y * (1.0/x12 - 1.0/x01) - np.roll(y, -1, axis=0) * (x01 / (x02*x12))
	
	d[0,...] = y[0,...] * (x01[1,...]+x02[1,...])/(x01[1,...]*x02[1,...]) - y[1,...] * x02[1,...]/(x01[1,...]*x12[1,...]) + y[2,...] * x01[1,...]/(x02[1,...]*x12[1,...])
	n2 = n-2
	d[n-1,...] = -y[n-3,...] * x12[n2,...]/(x01[n2,...]*x02[n2,...]) + y[n-2,...] * x02[n2,...]/(x01[n2,...]*x12[n2,...]) - y[n-1,...] * (x02[n2,...]+x12[n2,...]) / (x02[n2,...]*x12[n2,...])

	return d