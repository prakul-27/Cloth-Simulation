"""
This file is part of the submission "Homogenized Yarn-Level Cloth".
With this code, we provide our implementation of our algorithm for monotone
piecewise bicubic interpolation.

Inputs:
    tu          vector of length M of control points along u 
    tv          vector of length N of control points along v
    p           N x M array of control point values
    mu          N x M array of derivatives along u
    mv          N x M array of derivatives along v
    muv         N x M array of mixed derivatives
    extrapol    flag to indicate if first derivatives should be constrained
                for monotone extrapolation on the boundary

Usage:
    mu, mv, muv = monotone_bicubic(tu, tv, p, mu, mv, muv, extrapol)

    The u and v directions correspond to x and y in our supplementary
    explanation respectively.

    Note that the extrapol flag is only used to constrain the derivatives mu
    and mv on the boundaries. As discussed in our supplementary document, we
    additionally set muv = 0 on the boundary for extrapolation.

    We provide equation numbers from our supplementary document.
"""

import numpy as np

def constraint1(t, p):
    """Eq. (S36)"""
    return np.array([0, 3*(p[1]-p[0])/(t[1]-t[0])])

def enforce(c, m):
    """Project m to the range [c.min, c.max]"""
    if c is None:
        return m
    if m < c.min():
        return c.min()
    if m > c.max():
        return c.max()
    return m

def sweep1u(tu, p, mu):
    """Enforce (S36) on mu"""
    mu = mu.copy()
    for j in range(p.shape[0]):
        for i in range(p.shape[1]):
            if i > 0:
                mu[j,i] = enforce(constraint1(tu[i-1:i+1], p[j,i-1:i+1]), mu[j,i])
            if i+1 < p.shape[1]:
                mu[j,i] = enforce(constraint1(tu[i:i+2], p[j,i:i+2]), mu[j,i])
    return mu

def sweep1v(tv, p, mv):
    """Enforce (S36) on mv"""
    return sweep1u(tv, p.T, mv.T).T

def limit(d):
    """sign(d) * inf, such that it is 0 for d==0"""
    return np.sign(d)*np.inf if np.sign(d) != 0 else np.inf

def constraint2a(tv, p, i, j, extrapol):
    """Enforce (S40)"""
    d = p[j,i+1]-p[j,i]
    u = 3*d/(tv[j]-tv[j-1]) if j > 0 else (0 if extrapol else limit(d))
    l = -3*d/(tv[j+1]-tv[j]) if j+1 < len(tv) else (0 if extrapol else -limit(d))
    if np.isnan(u):
        print(i,j)
        print(u, d, np.sign(d)*np.inf)
    return np.array([l, u])

def absmax(a, b):
    """Return value farther from zero"""
    return a if abs(a) > abs(b) else b

def constraint2b(tu, tv, p, mu, i, j, extrapol):
    """Enforce (S41)"""
    d = p[j,i+1]-p[j,i]
    c = 3*d - (tu[i+1]-tu[i])*absmax(mu[j,i],mu[j,i+1])
    u = c/(tv[j]-tv[j-1]) if j > 0 else (0 if extrapol else limit(d))
    l = -c/(tv[j+1]-tv[j]) if j+1 < len(tv) else (0 if extrapol else -limit(d))
    return np.array([l, u])

def enforce2(c, m0, m1):
    """Find smallest change to (m0, m1) so that m1 - m0 lies in [c.min, c.max]  (Fig. S6)"""
    l = c.min()
    u = c.max()
    if m1 - m0 > u:
        if m0 + m1 > u:
            m0 = max(m0, 0)
            m1 = m0 + u
        elif m0 + m1 < -u:
            m1 = min(m1, 0)
            m0 = m1 - u
        else:
            a = (m0 + m1)/2
            m1 = a + u/2
            m0 = a - u/2
    if m1 - m0 < l:
        if m0 + m1 > -l:
            m1 = max(m1, 0)
            m0 = m1 - l
        elif m0 + m1 < l:
            m0 = min(m0, 0)
            m1 = m0 + l
        else:
            a = (m0 + m1)/2
            m1 = a + l/2
            m0 = a - l/2
    return m0, m1

def sweep2u(tu, tv, p, mu, mv, extrapol):
    """Enforce (S40), (S41) on mv"""
    mv = mv.copy()
    for j in range(p.shape[0]):
        for i in range(p.shape[1]-1):
            mv[j,i], mv[j,i+1] = enforce2(constraint2a(tv, p, i, j, extrapol=extrapol), mv[j,i], mv[j,i+1])
            mv[j,i], mv[j,i+1] = enforce2(constraint2b(tu, tv, p, mu, i, j, extrapol=extrapol), mv[j,i], mv[j,i+1])
        for i in reversed(range(p.shape[1]-1)):
            mv[j,i], mv[j,i+1] = enforce2(constraint2a(tv, p, i, j, extrapol=extrapol), mv[j,i], mv[j,i+1])
            mv[j,i], mv[j,i+1] = enforce2(constraint2b(tu, tv, p, mu, i, j, extrapol=extrapol), mv[j,i], mv[j,i+1])
    return mv

def sweep2v(tu, tv, p, mu, mv, extrapol):
    """Enforce (S40), (S41) on mu"""
    return sweep2u(tv, tu, p.T, mv.T, mu.T, extrapol=extrapol).T

def constraint3a(tu, tv, p, mu, mv, i, j):
    """Eq (S42a)"""
    if j+1 == p.shape[0] or i+1 == p.shape[1]:
        return None
    h = tu[i+1]-tu[i]
    k = tv[j+1]-tv[j]
    l = -3*mu[j,i]/k
    u = 3*((mv[j,i+1]-mv[j,i])/h + 3*(p[j,i+1]-p[j,i])/(h*k) - mu[j,i]/k)
    return np.array([l, u])

def constraint3b(tu, tv, p, mu, mv, i, j):
    """Eq (S42c)"""
    if j == 0 or i+1 == p.shape[1]:
        return None
    h = tu[i+1]-tu[i]
    k = tv[j]-tv[j-1]
    l = 3*((mv[j,i+1]-mv[j,i])/h - 3*(p[j,i+1]-p[j,i])/(h*k) + mu[j,i]/k)
    u = 3*mu[j,i]/k
    return np.array([l, u])

def constraint3c(tu, tv, p, mu, mv, i, j):
    """Eq (S42b)"""
    if j+1 == p.shape[0] or i == 0:
        return None
    h = tu[i]-tu[i-1]
    k = tv[j+1]-tv[j]
    l = -3*mu[j,i]/k
    u = 3*((mv[j,i]-mv[j,i-1])/h + 3*(p[j,i]-p[j,i-1])/(h*k) - mu[j,i]/k)
    return np.array([l, u])

def constraint3d(tu, tv, p, mu, mv, i, j):
    """Eq (S42d)"""
    if j == 0 or i == 0:
        return None
    h = tu[i]-tu[i-1]
    k = tv[j]-tv[j-1]
    l = 3*((mv[j,i]-mv[j,i-1])/h - 3*(p[j,i]-p[j,i-1])/(h*k) + mu[j,i]/k)
    u = 3*mu[j,i]/k
    return np.array([l, u])

def sweep3u(tu, tv, p, mu, mv, muv):
    """Enforce (S42a-d)"""
    muv = muv.copy()
    for j in range(p.shape[0]):
        for i in range(p.shape[1]):
            muv[j,i] = enforce(constraint3a(tu, tv, p, mu, mv, i, j), muv[j,i])
            muv[j,i] = enforce(constraint3b(tu, tv, p, mu, mv, i, j), muv[j,i])
            muv[j,i] = enforce(constraint3c(tu, tv, p, mu, mv, i, j), muv[j,i])
            muv[j,i] = enforce(constraint3d(tu, tv, p, mu, mv, i, j), muv[j,i])
    return muv

def sweep3v(tu, tv, p, mu, mv, muv):
    """Enforce (S42a-d)"""
    return sweep3u(tv, tu, p.T, mv.T, mu.T, muv.T).T

def monotone_bicubic(tu, tv, p, mu, mv, muv, extrapol=True):
    mu1 = sweep1u(tu, p, mu)
    mv1 = sweep1v(tv, p, mv)
    mv2 = sweep2u(tu, tv, p, mu1, mv1, extrapol=extrapol)
    mu2 = sweep2v(tu, tv, p, mu1, mv1, extrapol=extrapol)
    muv3 = sweep3u(tu, tv, p, mu2, mv2, muv)
    muv3 = sweep3v(tu, tv, p, mu2, mv2, muv)
    return mu2, mv2, muv3

