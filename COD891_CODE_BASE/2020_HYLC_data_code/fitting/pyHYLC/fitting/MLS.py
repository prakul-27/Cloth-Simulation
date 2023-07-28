import numpy as np

# bump function
def bbump(X, x0, kernel_size):
    d = (X - x0) / kernel_size
    vals = np.zeros_like(X)
    in_range = np.logical_and(d > -1 + 1e-10, d < 1 - 1e-10)
    vals[in_range] = np.exp(-1.0 / (1.0 - d[in_range]**2))
    return vals


# gaussian bump
def gbump(X, x0, kernel_size):
    d = (X - x0)
    vals = np.zeros_like(X)
    sig = kernel_size / 3.0
    in_range = np.logical_and(d > -3*sig, d < 3*sig)
    vals[in_range] = np.exp(-d[in_range]**2 / (2*sig**2)) / (np.sqrt(2*np.pi*sig**2))
    return vals

# evaluate MLS function given data X->Y at point x
def eval1D( X, Y, x, kernel_size, degree=3, W=None, bumpfun="bump"):
    assert(bumpfun in ["bump", "gaussian"])
    if bumpfun == "bump":
        w = bbump(X, x, kernel_size)
    else:
        w = gbump(X, x, kernel_size)
    in_range = w > 0
    V = np.vander(X[in_range], N=degree+1, increasing=True)

    w = w[in_range,np.newaxis]
    if not W is None:
        w *= W[in_range,np.newaxis]

    wV = w * V
    wVTV = np.einsum("ki,kj->ij",wV,V,optimize=True)
    wVTY = np.einsum("ki,k->i",wV,Y[in_range],optimize=True)
    c = np.linalg.solve(wVTV,wVTY)
    return x**np.array(range(0,degree+1)) @ c

# vandermonde matrix for 2d polynomials
def vander2d(X0, X1, N, only_mixed=False):
    coeffs = []
    for p in range(N):
        for q in range(N):
            if (p+q > N):
                continue
            if only_mixed and (p == 0 or q == 0):
                continue
            coeffs.append((p,q))

    V = np.zeros((len(X0),len(coeffs)))
    for i in range(len(coeffs)):
        (p,q) = coeffs[i]
        V[:,i] = X0**p * X1**q
    return coeffs, V

# evaluate MLS function given data (X0,X1)->Y at point (x0,x1)
def eval2D(X0, X1, Y, x0, x1, kernel_size, degree=3, W=None, bumpfun="bump"):
    assert(bumpfun in ["bump", "gaussian"])
    if bumpfun == "bump":
        w0 = bbump(X0, x0, kernel_size[0])
        w1 = bbump(X1, x1, kernel_size[1])
    else:
        w0 = gbump(X0, x0, kernel_size[0])
        w1 = gbump(X1, x1, kernel_size[1])
    w = w0 * w1
    in_range = w > 0

    w = w[in_range,np.newaxis]
    if len(w) == 0:
        print("WARNING KERNEL DOES NOT CONTAIN ANY POINTS")
        return 0

    if not W is None:
        w *= W[in_range,np.newaxis]
    
    coeffs, V = vander2d(X0[in_range], X1[in_range], N=degree+1, only_mixed=False)

    wV = w * V
    wVTV = np.einsum("ki,kj->ij",wV,V,optimize=True)
    wVTY = np.einsum("ki,k->i",wV,Y[in_range],optimize=True)
    c = np.linalg.solve(wVTV,wVTY)

    val = 0
    for i in range(len(coeffs)):
        (p,q) = coeffs[i]
        val += x0**p * x1**q * c[i]
    return val