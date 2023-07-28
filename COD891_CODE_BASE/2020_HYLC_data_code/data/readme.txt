Here, we provide the raw output from our microsolver.

Each pattern has a folder with an info-file and a data-file.

The "info.txt"-file contains the area of the periodic tile and the area density (computed as total mass / area).
The "data.txt"-file contains rows [s_x, s_a, s_y, II_00, II_01, II_11, Psi] for each experiment mapping input strains to output energy density Psi.

The order of experiments in the data file is as follows:
150x strains along s_x
150x strains along s_a
150x strains along s_y
150x strains along II_00
150x strains along II_11
(50x50)x strains along s_x s_a
(50x50)x strains along s_x s_y
(50x50)x strains along s_x II_00
(50x50)x strains along s_x II_11
(50x50)x strains along s_a s_y
(50x50)x strains along s_a II_00
(50x50)x strains along s_a II_11
(50x50)x strains along s_y II_00
(50x50)x strains along s_y II_11

Units:
s_x, s_a, s_y 		unitless
II_00, II_01, II_11 	1/m
Psi			N/m == kg/s^2
area			m^2
area density		kg/m^2	


The data is released under the MIT license (see LICENSE).
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
