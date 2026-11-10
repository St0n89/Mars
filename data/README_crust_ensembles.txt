This file collection contains models for the Martian crust derived from the inversion of Ps-receiver functions
of marsquakes S0173a, S0183a and S0235b as presented in Fig. 2 of 
Knapmeyer-Endrun, Panning, Bissig, Joshi, Khan et al.: Thickness and structure of the martian crust from InSight seismic data, Science, 373 (6553), 2021.
When using or displaying these models, please cite the above reference.

There are 8 data sets in total, containing the results of inversion methods A and B for low-frequency (LF) and high-frequency (HF)
receiver functions and using a 2-layer and a 3-layer crust, respectively. Each of these data sets contains 5000 models from the
distributions of models that can fit the data obtained in the inversions. The name of each data set contains information on the 
inversion method, the number of layers, and the type of receiver function data set used (LF or HF).

For inversion method A, a short description of the file contents is given at the start of each file. The models extend to about
400 km depth; at greater depth, the model TAYAK (Khan et al., PEPI, 258, 28-42, 2016, Smrekar et al., SSR, 215, 3, 2019), which
is also included in this distribution, was used. The models were provided by Felix Bissig and Amir Khan (ETH Zurich).

For inversion method B, the following parameters are provided for each layer: thickness (d), S-wave velocity (v), and vp/vs ratio
(r). The number following each parameter denotes the layer number, e.g. d2 denotes the thickness of layer 2. The models end with
a half-space for the mantle, which does not have a thickness parameter and is defined only by the S-wave velocity and vp/vs ratio.
Densities are not provided as they were not subject of the inversion, but calculated from the P-wave velocity via Birch's law
(Birch, JGR, 66, 2199-2224, 1961). The last column for each model denotes the L2 misfit value. The models were provided by
Rakshit Joshi (MPS Goettingen).
