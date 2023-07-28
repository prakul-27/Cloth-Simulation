import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches

bottom = cm.get_cmap('OrRd', 128)
top = cm.get_cmap('Blues_r', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
cmaperr = ListedColormap(newcolors)

symlog = lambda X: np.sign(X) * np.log(np.abs(X)+1)
symexp = lambda X: np.sign(X) * (np.exp(np.abs(X))-1)

strain_names = ["$s_x$","$s_a$","$s_y$","$II_{xx}$","$II_{yy}$"]

def plot1D(coeffs_funs, Xd, Yd, i, k, ax=None, n_tst=300, p_extend=0.1, ylabels=True, log=False):
    if ax is None:
        fig,ax = plt.subplots()

    squash = lambda Y: symlog(Y) if log else Y

    # plot ground truth data
    ax.plot(Xd[:,k],squash(Yd), "o", label='data')
        
    # get data 1D strain range and extend it
    rge = [np.min(Xd[:,k]),np.max(Xd[:,k])]
    rge = [rge[0] - p_extend * (rge[1]-rge[0]), rge[1] + p_extend * (rge[1]-rge[0])]
    
    # create, evaluate and plot test data
    C0,G,R,g,r,psi = coeffs_funs
    X_tst = np.zeros((n_tst, 6))
    X_tst[:,k] = np.linspace(rge[0],rge[1],n_tst)
    Y_tst = psi(X_tst)
    ax.plot(X_tst[:,k], squash(Y_tst), "C2-", label='fit')
    ax.plot(G[i][0],squash(C0 + G[i][1]), 'k.', label='control points')
        
    ax.set_xlabel(strain_names[i])
    # ax.set_title("$%s$" % k_to_str(i))
    if ylabels:
        ax.set_ylabel(r'log$(\Psi)$' if log else r"$\Psi$",rotation=0, labelpad=20)

def k_to_str(k,k1=-1):
    return "f_{%s%s}" % (['1','2','3','x','','y'][k], ['','1','2','3','x','','y'][k1+1])

def plot2D(coeffs_funs, Xd, Yd, i, coords, n_per_2D, figax=None, n_tst=15, p_extend=0.1, log=False,
                    cmap="Spectral_r", cmaperr=cmaperr, shading='flat'):
    if figax is None:
        fig, axs = plt.subplots(1, 5, figsize=(19.5, 4), sharex=True, sharey=True)
    else:
        fig, axs = figax
        
    k0, k1 = coords

    squash = lambda Y: symlog(Y) if log else Y
    
    vmin,vmax=None,None
    # Create a Rectangle patches
    a0,a1=Xd[:, k0].min(),Xd[:, k0].max()
    b0,b1=Xd[:, k1].min(),Xd[:, k1].max()
    rects = [patches.Rectangle((a0,b0),a1-a0,b1-b0,linewidth=1,edgecolor='0.5',facecolor='none') for n in range(3)]
    
    # plot fit
    C0,G,R,g,r,psi = coeffs_funs
    tu,tv = R[i][0],R[i][1]
    
    # get data 2D strain range and extend it
    rge0 = [np.min(Xd[:,k0]),np.max(Xd[:,k0])]
    rge0 = [rge0[0] - p_extend * (rge0[1]-rge0[0]), rge0[1] + p_extend * (rge0[1]-rge0[0])]
    rge1 = [np.min(Xd[:,k1]),np.max(Xd[:,k1])]
    rge1 = [rge1[0] - p_extend * (rge1[1]-rge1[0]), rge1[1] + p_extend * (rge1[1]-rge1[0])]

    # create, evaluate and plot test data
    X0 = np.linspace(rge0[0],rge0[1],n_tst)
    X1 = np.linspace(rge1[0],rge1[1],n_tst)
    X0 -= X0[np.argmin(np.abs(X0))]
    X1 -= X1[np.argmin(np.abs(X1))]
    X0,X1 = np.meshgrid(X0,X1)
    
    X_tst = np.array([[0.0,0,0,0,0,0]]*n_tst*n_tst)
    X_tst[:,k0] = X0.reshape(-1)
    X_tst[:,k1] = X1.reshape(-1)
    
    Y_tst = psi(X_tst)

    vmin = squash(Y_tst.min())
    vmax = squash(Y_tst.max())

    Y_tst = Y_tst.reshape(n_tst,n_tst)
    im = axs[3].pcolormesh(X0,X1, squash(Y_tst), vmin=vmin,vmax=vmax,cmap=cmap,shading=shading)
    fig.colorbar(im, ax=axs[3])
    axs[3].set_title(r"$\Psi_{fit}$")


    # residual
    i0 = [0,1,2,3,4,4][k0]
    i1 = [0,1,2,3,4,4][k1]
    Yr = squash(g[i0](X_tst).reshape(n_tst,n_tst)+g[i1](X_tst).reshape(n_tst,n_tst))
    im = axs[1].pcolormesh(X0,X1, Yr, vmin=vmin,vmax=vmax,cmap=cmap,shading=shading)
    fig.colorbar(im, ax=axs[1])
    axs[1].set_title("$%s+%s$"%(k_to_str(i0),k_to_str(i1)))

    Yr = squash(r[i](X_tst).reshape(n_tst,n_tst))
    vmaxr = np.abs(Yr).max()
    im = axs[2].pcolormesh(X0,X1, Yr, vmin=-vmaxr,vmax=vmaxr,cmap=cmaperr,shading=shading)
    fig.colorbar(im, ax=axs[2])
    axs[2].set_title("$%s$"%k_to_str(i0,i1))

    # Add the patch to the Axes
    for i in range(len(rects)):
        axs[i+1].add_patch(rects[i])
        
    # plot data (after for vmin vmax)
    XX = Xd[:, k0].reshape(n_per_2D,n_per_2D)
    YY = Xd[:, k1].reshape(n_per_2D,n_per_2D)
    im = axs[0].pcolormesh(XX, YY, squash(Yd).reshape(n_per_2D,n_per_2D),
                          vmin=vmin,vmax=vmax,cmap=cmap,shading=shading)
    fig.colorbar(im, ax=axs[0])
    axs[0].set_title("data")

    # plot relative error
    # reevaluate Ydtst from Xd, compute err = Ydtest - Yd/(abs(Yd) + 1e-5)
    Yd_tst = psi(Xd)
    # err = (Yd_tst - Yd)/(np.abs(Yd) + 1e-5)
    # axs[4].set_title("Rel.Err.")
    err = squash(Yd_tst - Yd)
    axs[4].set_title(r"$\Psi_{fit}$ - data")
    vmaxr = np.abs(err).max() * 0.8
    im = axs[4].pcolormesh(XX,YY, err.reshape(n_per_2D,n_per_2D), vmin=-vmaxr,vmax=vmaxr,cmap=cmaperr,shading=shading)
    fig.colorbar(im, ax=axs[4])

    # control point grid
    for ax in axs:
        for u in tu:
            ax.axvline(u, c='k', linestyle="-" if u == 0 else "--", alpha=0.2 if u == 0 else 0.125, linewidth=1)
        for v in tv:
            ax.axhline(v, c='k', linestyle="-" if v == 0 else "--", alpha=0.2 if v == 0 else 0.125, linewidth=1)

        TU,TV = np.meshgrid(tu,tv)
        ax.plot(TU,TV,'k.',alpha=0.2)
                
        ax.set_xlim(rge0)
        ax.set_ylim(rge1)
        
    for i,ax in enumerate(axs):
        ax.set_xlabel(strain_names[i0])
        ax.set_ylabel(strain_names[i1])
    rhs = r"symlog(\Psi)" if log else r"\Psi"
    if figax is None:
        fig.suptitle(r"%s, %s $\rightarrow %s$" % (strain_names[i0],strain_names[i1],rhs))
    else:
        pass
        axs[0].annotate(r"Strains: %s, %s" % (strain_names[i0],strain_names[i1]), xy=(0.0, 1.2), 
                xycoords='axes fraction',
                size='large', ha='center', va='top')
        # plt.annotate(r"%s, %s $\rightarrow %s$" % (strain_names[i0],strain_names[i1],rhs), xy=(1+0.5, 1.2), 
        #         xycoords='axes fraction',
        #         size='large', ha='center', va='top')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
