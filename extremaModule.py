def extrema(x, max = True, min = True, strict = False, withend = False):
    """
    This function will index the extrema of a given array x.
	
    Options:
		max		If true, will index maxima
		min		If true, will index minima
		strict		If true, will not index changes to zero gradient
		withend	If true, always include x[0] and x[-1]
	
    This function will return a tuple of extrema indicies and values
    """
	
    # This is the gradient
    from numpy import zeros
    dx = zeros(len(x))
    from numpy import diff
    dx[1:] = diff(x)
    dx[0] = dx[1]
	
    # Clean up the gradient in order to pick out any change of sign
    from numpy import sign
    dx = sign(dx)
	
    # define the threshold for whether to pick out changes to zero gradient
    threshold = 0
    if strict:
    	threshold = 1
		
    # Second order diff to pick out the spikes
    d2x = diff(dx)
	
    if max and min:
     	d2x = abs(d2x)
    elif max:
 	d2x = -d2x
	
    # Take care of the two ends
    if withend:
	d2x[0] = 2
 	d2x[-1] = 2
	
    # Sift out the list of extremas
    from numpy import nonzero
    ind = nonzero(d2x > threshold)[0]
	
    return ind, x[ind]
