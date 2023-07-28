import numpy as np
import os
import json

def read(folder):
    data = np.loadtxt(os.path.join(folder,"data.txt"))
    X = data[:,:6] # [sx, sa, sy, II00, II01=0, II11]
    Y = data[:,6]  # Psi
    
    normX_scale = np.max(np.abs(X), axis=0)
    normX_scale[4] = 1 # no normalization needed for II_01 = 0

    X = X / normX_scale

    info = np.loadtxt(os.path.join(folder,"info.txt"))
    area_density = info[1]

    return X, Y, normX_scale, area_density


def material_dict(filename, area_density, ixs_2D, coeffs, normX_scale):#, strain_range):
    mat = {}
    mat["density"] = area_density
    mat["strain scale"] = normX_scale.tolist()
    mat["coeffs"] = {}
    mat["coeffs"]["const"] = coeffs[0]
    mat["coeffs"]["1D"] = []
    mat["coeffs"]["2D"] = []
    
    # 1d
    for k in range(5):
        (t,p,m) = coeffs[1][k]
        coeffjson = {}
        coeffjson["k"] = k
        coeffjson["t"] = t.tolist()
        coeffjson["p"] = p.tolist()
        coeffjson["m"] = m.tolist()
        mat["coeffs"]["1D"].append(coeffjson)

    # 2d
    for k in range(min(len(ixs_2D), len(coeffs[2]))):
        k0, k1 = ixs_2D[k]
        (tu, tv, p,mu,mv,muv) = coeffs[2][k]
        coeffjson = {}
        coeffjson["k0"] = int(k0)
        coeffjson["k1"] = int(k1)
        coeffjson["tu"] = tu.tolist()
        coeffjson["tv"] = tv.tolist()
        coeffjson["p"] = p.tolist()
        coeffjson["mu"] = mu.tolist()
        coeffjson["mv"] = mv.tolist()
        coeffjson["muv"] = muv.tolist()
        mat["coeffs"]["2D"].append(coeffjson)

    if not filename is None:
        with open(filename, 'w') as thefile:
            json.dump(mat, thefile, sort_keys=False, indent=2)
    else:
        return mat

