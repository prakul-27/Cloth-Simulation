Here, we provide python code for fitting the raw strain-energy-pair data.
The code is located in the 'pyHYLC' directory, and we provide a jupyter notebook
'fitting.ipynb' for loading and fitting data, plotting the results, and saving
the coefficients in a json file.

A small explanation about the json format:
Each material is given as a json file with the following entries:
"density":      is the area density of the material in kg/m^2
"strain shift": can be ignored. (refers to the fact that sx and sy are defined to be 0-centered)
"coeffs":
   "const": the constant part Psi0
   "1D": lists of k,t,p,m where k denotes strain index (starting at 0 with the order: sx,sa,sy,IIxx,IIyy), and t,p,m are the spline control point locations, values, and tangents respectively.
   "2D": lists of k0,k1,tu,tv,p,mu,mv,muv. k0,k1 denote the strain indices like for 1D and tu,tv,p,mu,mv,muv are 2D spline control points locations (tu,tv), values, first (mu,mv) and mixed derivatives (muv) in the u and v directions respectively. The arrays p,mu,mv,muv are flattened row-major, i.e. the coefficents corresponding to the node with index 2 along u and index 3 along v are given by tu[2],tv[3],p[3*len(tu)+2],etc.
Thus, the 1D spline with entry k=0 corresponds to f_1 in the paper, k=3 and k=4 correspond to f_x and f_y respectively, and similarly for 2D splines.
We refer to the supplementary document Equation (S26a-c) for the verbose definition of our energy model.

In the code we have (partially) adopted a different notation of spline coefficients. For completeness, the following equivalences hold: 
1D: t_i = x_i, p_i = p_i, m_i = p^x_i
2D: {t_u}_ij = x_ij, {t_v}_ij = y_ij, p_ij = p_ij, mu_ij = p^x_ij, mv_ij = p^y_ij, muv_ij = p^xy_ij


The code is released under the MIT license (see LICENSE).
If you use our code or data, please consider citing us:
@article{sperl2020hylc,
  author    = {Sperl, Georg and Narain, Rahul and Wojtan, Chris}
  title     = {Homogenized Yarn-Level Cloth},
  journal   = {ACM Transactions on Graphics (TOG)},
  number    = {4},
  volume    = {39},
  year      = {2020},
  publisher = {ACM}
}
