import numpy as np

def monotone_cubic(spl):
    (t,p,m) = [np.copy(arr) for arr in spl]
    # https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
    # NOTE wikipedia notation is  1-indexed

    dk = (p[1:]-p[:-1])/(t[1:]-t[:-1])

    # prohibit change of tangent sign between neighbors
    for k in range(1, len(t)-1):
        if dk[k] * dk[k-1] < 0:
            m[k] = 0

    for k in range(0, len(t)-1):
        # enforce monotonicity between flat neighbours
        if abs(p[k+1]-p[k]) < 1e-8: # APPROX
            m[k] = 0
            m[k+1] = 0
            continue

        # alpha = m[:-1] / dk
        # beta = m[1:] / dk
        alpha = m[k] / dk[k]
        beta = m[k+1] / dk[k]

        # if the secant sign does not match with the tangent signs
        if alpha < 0:
            m[k] = 0
        if beta < 0:
            m[k+1] = 0

        # prevent overshoot METHOD 1
        a2b2 = alpha**2+beta**2
        if a2b2 > 9:
            tau = 3/np.sqrt(a2b2)
            m[k] = tau * alpha * dk[k]
            m[k+1] = tau * beta * dk[k]
            
    return (t,p,m)
    
# evaluate Cubic Hermite Spline given the (t,p,m) data
def evalCHS(spl, X, ext=0):
    # unpack spline
    (t,p,m) = spl
 
    # NOTE could write in matrix form for array evaluation
    # X = X * 0 + -0.4
    # get interval per x
    Y = np.empty_like(X)
    for k,x in enumerate(X):
        # if k >= 1:
        #     Y[k:] = 0
        #     break
        i = np.searchsorted(t, x, side="right")
        # print(x, i)
        if i == 0 or i == len(t):
            if i == 0:
                tt, pp, mm = t[0], p[0], m[0]#, t[1]-t[0]
            else:
                tt, pp, mm = t[-1], p[-1], m[-1]#, t[-1]-t[-2]
            if ext == 1:
                mm = 0 # constant
            elif ext == 2:
                pp, mm = 0, 0
            Y[k] = pp + mm * (x - tt) 
        else:
            # p(t) = (2t^3 - 3t^2 + 1) p0 + (t^3 - 2t^2 + t) m0 + (-2t^3 + 3t^2)p1 + (t^3 - t^2) m1
            # dp(t) = (6t^2 - 6t^1 ) p0 + (3t^2 - 4t + 1) m0 + (-6t^2 + 6t)p1 + (3t^2 - 2t) m1
            # ddp(t) = (12 t - 6 ) p0 + (6t - 4) m0 + (-12t + 6)p1 + (6t - 2) m1
            t1, t0 = t[i], t[i-1]
            p1, p0 = p[i], p[i-1]
            m1, m0 = m[i], m[i-1]
 
            tt = (x - t0) / (t1 - t0) # pos in 0,1 interval
            h00 = (2 *( tt**3 )- 3* (tt**2) + 1)
            h10 = (tt**3 - 2* (tt**2) + tt)
            h01 = (-2* (tt**3) + 3*(tt**2))
            h11 = ((tt**3) - (tt**2))
            Y[k] = (h00 * p0
                + h10 * (t1 - t0) * m0
                + h01 * p1
                + h11 * (t1 - t0) * m1)
 
    return Y

# evaluate 2D Cubic Hermite Spline given the (tu, tv, p, mu, mv, muv) data
def evalCHS2D(spl, X0, X1, ext=0, dx=0, dy=0, order=3):
    # unpack spline
    (tu, tv, p, mu, mv, muv) = spl
 
    # NOTE could write in np matrix form for array evaluation maybe
 
    nu = len(tu)
    # nv = len(tv)
 
    if not order in [0,1,3]:
        print("no kernel for order", order," - defaulting to 3")
        order = 3
    # print("ORDER",order, len(tu),len(tv))
    if order == 3:
        M = np.array([
            [2,-2,1,1],
            [-3,3,-2,-1],
            [0,0,1,0],
            [1,0,0,0]
        ])
    elif order == 1:
        M = np.array([
            [0,0,0,0],
            [0,0,0,0],
            [-1,1,0,0],
            [1,0,0,0]
        ])
    elif order == 0:
        M = np.array([
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0.5,0.5,0,0]
        ])

    #Y = np.empty_like(X0, dtype=np.float)
    #--change by me---
    Y = np.empty_like(X0, dtype=float)
    for k,(x0,x1) in enumerate(zip(X0,X1)):
        j = np.searchsorted(tu, x0, side="right")
        i = np.searchsorted(tv, x1, side="right")

        extrapol_u = False
        extrapol_u_right = False
        extrapol_v = False
        extrapol_v_right = False
        if i <= 0 or i >= len(tv):
            extrapol_v = True
            extrapol_v_right = i >= len(tv)
            i = max(1,min(len(tv)-1,i))

        if j <= 0 or j >= len(tu):
            extrapol_u = True
            extrapol_u_right = j >= len(tu)
            j = max(1,min(len(tu)-1,j))

        deltau = (tu[j] - tu[j-1])
        if extrapol_u_right:
            u = (x0 - tu[j]) / deltau
        else:
            u = (x0 - tu[j-1]) / deltau
        deltav = (tv[i] - tv[i-1])
        if extrapol_v_right:
            v = (x1 - tv[i]) / deltav
        else:
            v = (x1 - tv[i-1]) / deltav
         
        i = i - 1 # shift indexing to left
        j = j - 1

        # order BL=00,BR=01,TL=10,TR=11
        p00, p01, p10, p11 = p[i*nu + j], p[i*nu + j+1], p[(i+1)*nu + j], p[(i+1)*nu + j+1]
        mu00, mu01, mu10, mu11 = mu[i*nu + j], mu[i*nu + j+1], mu[(i+1)*nu + j], mu[(i+1)*nu + j+1]
        mv00, mv01, mv10, mv11 = mv[i*nu + j], mv[i*nu + j+1], mv[(i+1)*nu + j], mv[(i+1)*nu + j+1]
        muv00, muv01, muv10, muv11 = muv[i*nu + j], muv[i*nu + j+1], muv[(i+1)*nu + j], muv[(i+1)*nu + j+1]
 

        # ext: 0 -> linear, 1 -> const, 2 -> 0

        deltas = 1.0
        if dx == 0:
            U = u**np.array([3,2,1,0])
        elif dx == 1:
            U = np.array([3 * u * u, 2 * u, 1.0, 0.0])
            deltas *= 1.0 / deltau
        elif dx == 2:
            U = np.array([6 * u, 2.0, 0.0 ,0.0])
            deltas *= 1.0 / deltau**2
        else:
            assert(False)

        if dy == 0:
            V = v**np.array([3,2,1,0])
        elif dy == 1:
            V = np.array([3 * v * v, 2 * v, 1.0, 0.0])
            deltas *= 1.0 / deltav
        elif dy == 2:
            V = np.array([6 * v, 2.0, 0.0 ,0.0])
            deltas *= 1.0 / deltav**2
        else:
            assert(False)

        Mu = np.copy(M)
        Mv = np.copy(M)
        lin = 1 if ext == 0 else 0  #if 1 then linear else constant
        if extrapol_u:
            if ext==2:
                Y[k] = 0
                continue
                
            if not extrapol_u_right: # if left side
                Mu = np.array([ # NOTE to self: M such that [U..] . M . [p0,p1,m0,m1]^T
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,lin,0], # p0 + u * m0/du
                    [1,0,0,0]
                ])
            else: # if right side
                Mu = np.array([
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,0,lin], # p1 + u * (m1)/du
                    [0,1,0,0]
                    # [0,1,0,-lin]
                ])
        if extrapol_v:
            if ext==2:
                Y[k] = 0
                continue

            if not extrapol_v_right:
                Mv = np.array([
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,lin,0],
                    [1,0,0,0]
                ])
            else:
                Mv = np.array([
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,0,lin],
                    [0,1,0,0]
                ])

        B = np.array([
            [p00,p10, mv00, mv10],
            [p01,p11, mv01, mv11],
            [mu00,mu10, muv00, muv10],
            [mu01,mu11, muv01, muv11]
        ])
        B = B * np.array([
            [1,1,deltav,deltav],
            [1,1,deltav,deltav],
            [deltau,deltau,deltau*deltav,deltau*deltav],
            [deltau,deltau,deltau*deltav,deltau*deltav],
        ])
 
        # U.T M B M.T V
        Y[k] = deltas * np.einsum("i,ij,jk,lk,l->",U,Mu,B,Mv,V)
 
    return Y