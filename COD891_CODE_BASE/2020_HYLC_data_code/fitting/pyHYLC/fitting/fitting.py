import numpy as np
from . import splines, MLS, monotone_bicubic

# vec = lambda *x: np.array(x,dtype=np.float)
#---change by me
vec = lambda *x: np.array(x,dtype=float)
symlog = lambda X: np.sign(X) * np.log(np.abs(X)+1)
symexp = lambda X: np.sign(X) * (np.exp(np.abs(X))-1)

def quasiconvexify1D(t, p, m, i_min, eps=1e-8):
    """
    march outward from minimum to ensure minimum increase in values and consistent tangents
    updates p and m (depending on which side of minimum)
    - p_i+1 >= p_i + (t_i+1 - t_i) * eps 
    - m_i >= eps
    """
    p = np.copy(p)

    # march outward left and right and project values p
    for i in range(i_min + 1, len(t)): # right of min
        d = np.abs(t[i]-t[i-1])
        p[i] = max(p[i], p[i-1] + d * eps) # B >= A + d eps
    for i in range(i_min - 1, -1, -1): # left of min
        d = np.abs(t[i]-t[i+1])
        p[i] = max(p[i], p[i+1] + d * eps)
    
    # project derivatives
    m[i_min] = 0
    for i in range(i_min + 1, len(t)): # right of min
        m[i] = max(m[i], eps)
    for i in range(i_min - 1, -1, -1): # left of min
        m[i] = min(m[i], -eps)
    return p, m


def fit_1D(X, Y, k, reg=0, symmetries=[]):
    """fit 1D splines for each 1D range"""

    # set spline ctrl point positions as linspaced within 95% of data min/max
    n = 11
    s = 0.95
    t = np.linspace(X[:,k].min()*s, X[:,k].max()*s, n)
    t -= t[np.argmin(np.abs(t))] # ensure passing through origin
        
    # smooth data using log-MLS
    logY = symlog(Y)
    kernel_size = ((X[:,k].max()-X[:,k].min())*0.15) # MLS kernel as 15% of span
    if k in symmetries:
        mls_base = lambda u: MLS.eval1D(X[:,k], logY, u, kernel_size,
                                    degree=2, bumpfun="gaussian") 
        mls = lambda u: 0.5*(mls_base(u)+mls_base(-u))
    else:
        mls = lambda u: MLS.eval1D(X[:,k], logY, u, kernel_size,
                                    degree=2, bumpfun="gaussian")
    f = lambda u: symexp(mls(u)) # final smoothed function
    
    # get values p and tangents m by evaluating and finite-differencing MLS
    p = np.empty_like(t) # spline ctrl point values
    m = np.empty_like(t) # spline ctrl point derivatives
    eps = np.min(t[1:]-t[:-1]) * 0.1
    for i in range(len(p)):
        p[i] = f(t[i]) # value
        m[i] = (f(t[i]+0.5*eps) - f(t[i]-0.5*eps)) / eps # derivative

    # quasiconvexify the ctrl points, marching outward from minimum i_min
    i0 = np.argmin(np.abs(t))
    if k < 3: # in-plane minimum assumed at 0
        i_min = i0
    else: # bending minimum found at minimum p
        i_min = np.argmin(p)
    p, m = quasiconvexify1D(t, p, m, i_min=i_min, eps=reg)

    # apply Fritsch-Carlson method for cubic monotone spline interpolation
    t,p,m = splines.monotone_cubic((t,p,m))

    # increasing extrapolation
    eps = 1e-8
    m[0] = min(m[0],-eps)
    m[-1] = max(m[-1],eps)

    # residualize by subtracting value at t=0 (because the constant C = Psi(0,...) is already fit)
    p -= p[i0]

    return (t,p,m)


def quasiconvexify2D(tu,tv,p, mu=None,mv=None,i_min=None,j_min=None, eps=1e-8):
    """
    march outward from minimum to ensure minimum increase in values and consistent tangents
    updates p and m (depending on which side of minimum)
    - p_i+1,j >= p_i,j + (tv_i+1 - tv_i) * eps 
    - p_i,j+1 >= p_i,j + (tu_j+1 - tu_j) * eps 
    - mu_i,j >= eps
    - mv_i,j >= eps
    """
    p = np.copy(p)

    # get axes
    i0 = np.argmin(np.abs(tv))
    j0 = np.argmin(np.abs(tu))
    assert(np.abs(tu[j0])<1e-6 and np.abs(tv[i0])<1e-6)

    if i_min is None:
        i_min = np.argmin(np.abs(tv))
    if j_min is None:
        j_min = np.argmin(np.abs(tu))

    # get ctrl point ixs sorted by distance to min, ensuring correct L1 traversal
    II,JJ = np.meshgrid(np.arange(len(tv)),np.arange(len(tu)))
    IJ = np.array(list(zip(II.reshape(-1),JJ.reshape(-1))))
    IJ = IJ[np.argsort([(i-i_min)**2 + (j-j_min)**2 for (i,j) in IJ])]

    # project each ctrl point from both sides
    for (i,j) in IJ:
        if (i,j) == (i_min,j_min):
            continue
        
        # assume p on axes is already set to correct value (ie 1D stuff)
        # so never overwrite but use the hopefully qconvex values
        if (i == i0 or j == j0):
            continue
        
        diru = np.sign(j-j_min)
        dirv = np.sign(i-i_min)
        B = p[i*len(tu)+j]
        A = p[(i-dirv)*len(tu) + j]
        d = np.linalg.norm(vec(tu[j],tv[i])-vec(tu[j],tv[i-dirv]))
        Bproj = max(B, A + d * eps)

        A = p[i*len(tu) + (j-diru)]
        d = np.linalg.norm(vec(tu[j],tv[i])-vec(tu[j-diru],tv[i]))
        Bproj = max(Bproj, A + d * eps)

        p[i*len(tu)+j] = Bproj

        if not mu is None and not mv is None:
            ix = i*len(tu)+j
            mu[ix] = max(mu[ix], eps) if diru > 0 else min(mu[ix], - eps)
            mv[ix] = max(mv[ix], eps) if dirv > 0 else min(mv[ix], - eps)
    if not mu is None and not mv is None:
        return p,mu,mv
    return p


def fit_2D(X, Y, k0, k1, reg=0, symmetries=[], stockinette=False):
    def resample(t0,t1,sym):
        # resample the interval [t0,0,t1] as linspaced [t0,0] and linspaced [0,t1]
        # NOTE: nleft, nright including 0
        # sym is False only for s_x and s_y, where we need less samples in compression

        s = 0.9
        if sym:
            nleft = 5
            nright = 5
        else:
            nleft = 4 if [k0,k1] == [0,2] else 2 # poisson gets 3 left of axes, others only need 1 because 0 residual 
            nright = 7
        
        if stockinette and k1 > 2: # stockinette bending terms
            s = 0.3
            if sym:
                nleft, nright = 3, 3
            else:
                nleft, nright = 2, 4

        t = np.concatenate([
            np.linspace(t0*s,0,nleft),
            np.linspace(0,t1*s,nright)[1:]
        ])
        return t

    # space control points tu,tv according to resampling function above
    tu = resample(X[:,k0].min(), X[:,k0].max(), sym=k0 not in [0,2])
    tv = resample(X[:,k1].min(), X[:,k1].max(), sym=k1 not in [0,2])

    # MLS
    # f(u,v) = symexp(mls(X0, X1, logY, u, v))
    kernel_size = [(X[:,k].max()-X[:,k].min())*0.15 for k in [k0,k1]]
    Ylog = symlog(Y)
    mls_base = lambda u,v: MLS.eval2D(X[:,k0], X[:,k1], Ylog, u, v,
        kernel_size, degree=2, W=None, bumpfun="gaussian")
    # symmetrize
    if k0 in symmetries and k1 in symmetries:
        mls = lambda u,v: 0.25*(mls_base(u,v) + mls_base(-u,v)
                                + mls_base(u,-v) + mls_base(-u,-v))
    elif k0 in symmetries:
        mls = lambda u,v: 0.5*(mls_base(u,v) + mls_base(-u,v))
    elif k1 in symmetries:
        mls = lambda u,v: 0.5*(mls_base(u,v) + mls_base(u,-v))
    else:
        mls = mls_base
    f = lambda u,v: symexp(mls(u,v)) # final interpolation function
    
    # set values p = f(tu, tv)
    # and finite-difference derivatives mu, mv, muv
    p = np.empty((len(tu)*len(tv))) 
    mu = np.empty_like(p)
    mv = np.empty_like(p)
    muv = np.zeros_like(p)
    eps= min(np.min(tu[1:]-tu[:-1]),np.min(tv[1:]-tv[:-1])) * 0.1 # FD-epsilon
    for i in range(len(tv)):
        for j in range(len(tu)):
            ix = i*len(tu) + j
            u, v = tu[j], tv[i]

            # value
            p[ix] = f(u,v)
            # derivative mu, mv
            mu[ix] = (f(u+0.5*eps,v) - f(u-0.5*eps,v)) / eps
            mv[ix] = (f(u,v+0.5*eps) - f(u,v-0.5*eps)) / eps
            # derivative muv
            muv[ix] = ((f(u+0.5*eps,v+0.5*eps) - f(u-0.5*eps,v+0.5*eps)) - (f(u+0.5*eps,v-0.5*eps) - f(u-0.5*eps,v-0.5*eps))) / (eps**2)


    # quasiconvexify by L1-marching (except poisson)
    if not (k0 == 0 and k1 == 2):
        # define minima
        minimum_u = 0 if k0 < 3 else None
        minimum_v = 0 if k1 < 3 else None
        min_type = int(minimum_u is None) + int(minimum_v is None)
        if min_type == 2: # undefined min, define via smallest knot value (NOTE: does not happen)
            ix_min = np.argmin(p)
            i_min = ix_min // len(tu)
            j_min = ix_min % len(tu)
        elif min_type == 1: # one undefined min, find as constrained smallest knot (NOTE: mixed in-plane/bending)
            if minimum_u is None:
                i_min = np.argmin(np.abs(tv))
                pp = np.array([p[i_min*len(tu)+j] for j in range(len(tu))]) # v slice
                j_min = np.argmin(pp)
            else:
                j_min = np.argmin(np.abs(tu))
                pp = np.array([p[i*len(tu)+j_min] for i in range(len(tv))]) # u slice
                i_min = np.argmin(pp)
        else: # defined min assumed as (0,0) (NOTE: in-plane deformation only)
            i_min = np.argmin(np.abs(tv))
            j_min = np.argmin(np.abs(tu))
            
        p,mu,mv = quasiconvexify2D(tu,tv,p,mu,mv,
            i_min=i_min, j_min=j_min,
            eps=reg)

    # Enforce increasing extrapolation
    eps = 1e-8
    mu, mv = mu.reshape(len(tv), len(tu)), mv.reshape(len(tv), len(tu))
    mu[:,0] = np.minimum(mu[:,0], -eps)
    mu[:,-1] = np.maximum(mu[:,-1], eps)
    mv[0,:] = np.minimum(mv[0,:], -eps)
    mv[-1,:] = np.maximum(mv[-1,:], eps)
    mu, mv = mu.reshape(-1), mv.reshape(-1)

    # Apply bicubic monotone algorithm
    [p,mu,mv,muv] = [arr.reshape(tv.size, tu.size) for arr in [p,mu,mv,muv]]
    mu, mv, muv = monotone_bicubic.monotone_bicubic(tu, tv, p, mu, mv, muv, extrapol=True)
    muv[[0,-1],:] = 0
    muv[:,[0,-1]] = 0
    [p,mu,mv,muv] = [arr.reshape(-1) for arr in [p,mu,mv,muv]]

    # Convert to residual: r = f(x,y) - f(x,0) - f(0,y) + f(0,0)
    # for p: subtract axes from all u and v slices and add origin
    # for mu: subtract axis from all u slices
    # for mv: subtract axis from all v slices
    # for muv: == muv

    # indices of origin/axes
    i0 = np.argmin(np.abs(tv))
    j0 = np.argmin(np.abs(tu))

    for i in range(len(tv)):
        if i == i0:
            continue
        for j in range(len(tu)):
            if j==j0:
                continue
            p[i*len(tu)+j] = (p[i*len(tu)+j] - p[i0*len(tu)+j]
                                    - p[i*len(tu)+j0] + p[i0*len(tu)+j0])
            mu[i*len(tu)+j] = mu[i*len(tu)+j] - mu[i0*len(tu)+j]
            mv[i*len(tu)+j] = mv[i*len(tu)+j] - mv[i*len(tu)+j0]

    for i in range(len(tv)):
        p[i*len(tu)+j0] = 0
        mu[i*len(tu)+j0] = 0
        mv[i*len(tu)+j0] = 0
    for j in range(len(tu)):
        p[i0*len(tu)+j] = 0
        mu[i0*len(tu)+j] = 0
        mv[i0*len(tu)+j] = 0

    # Except for the poisson mode (s_x, s_y) enforce 0 compression residual (due to noise in the data)
    if not (k0 == 0 and k1 == 2):
        compr_ixs = np.zeros((len(tv),len(tu))).astype(bool)
        if k0 in [0,2]:
            compr_ixs[:, tu < 0] = True
        if k1 in [0,2]:
            compr_ixs[tv < 0, :] = True
        compr_ixs = compr_ixs.reshape(-1)

        p[compr_ixs] = 0
        mu[compr_ixs] = 0
        mv[compr_ixs] = 0
        muv[compr_ixs] = 0
    
    spl = (tu, tv, p, mu, mv, muv)

    return spl


def fit(X, Y, ixs_1D, blocks_1D, ixs_2D, blocks_2D, stockinette=False):
    """
    Fit 1D and 2D splines from normalized strains X and energy densities Y.
    1D and 2D subdata are given through coordinates/slices ixs_1D/block_1D and ixs_2D/blocks_2D.
    """

    reg = 1e-3 # minimum absolute value of increase of quasiconvexity
    if stockinette:
        reg = 1e-2 # higher regularization for the stockinette

    symmetries = [1] # strain coords. around which to symmetrize MLS (i.e. coord of s_a)

    # fit constant energy minimum == Psi(0,0,0,0,0,0)
    ix = np.argmin(np.linalg.norm(X - np.array([0,0,0,0,0,0]),axis=1))
    C0 = Y[ix]
    
    # fit 1D terms
    coeffs_1d = []
    for k, block in enumerate(blocks_1D):
        spl = fit_1D(X[block], Y[block], ixs_1D[k], reg=reg, symmetries=symmetries)
        coeffs_1d.append(spl)
    
    # fit 2D terms
    coeffs_2d = []
    for k, block in  enumerate(blocks_2D):
        k0, k1 = ixs_2D[k]
        spl = fit_2D(X[block], Y[block], k0, k1, reg=reg, symmetries=symmetries, stockinette=stockinette)
        coeffs_2d.append(spl)
    
    coeffs = [C0, coeffs_1d, coeffs_2d]
    return coeffs


def coeffs_to_splines(coeffs, ixs_1D, ixs_2D):
    C0 = coeffs[0]
    G = coeffs[1] # 1D spline coeffs: G[k] = t,p,m
    R = coeffs[2] # 2D spline coeffs: R[i] = tu,tv,p,mu,mv,muv

    # individual spline evaluation functions
    g = [lambda Xnew,i=i,k=k: splines.evalCHS(coeffs[1][i], Xnew[:,k]) for i,k in enumerate(ixs_1D)]
    r = [lambda Xnew,i=i,k0=k0,k1=k1:
         splines.evalCHS2D(coeffs[2][i], Xnew[:,k0], Xnew[:,k1])
         for i,(k0,k1) in enumerate(ixs_2D)]
    
    # NOTE: this psi is a simplified version of the full model used for plotting
    # The total model simplifies to the sum of individual splines under the assumption that we only query pre-normalized strains Xnew with only single curvature as in the data (II00, 0, 0) or (0, 0, II11)
    # i.e. there is no need to normalize the input and compute principal curvatures here
    psi_simplified = lambda Xnew: (C0
                        + sum([g[i](Xnew) for i in range(len(g))])
                        + sum([r[i](Xnew) for i in range(len(r))])
                       )
    return C0,G,R,g,r,psi_simplified
