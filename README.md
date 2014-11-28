<<<<<<< HEAD
﻿TRT-SInterp
=======
TRT_SInterp
>>>>>>> origin/master
===========

Stochastic Interpretation of Thermal Response Test

## Description
TRT-SInterp is a Matlab code designed to interpret a thermal response test in a deterministic or stochastic framework. The program treats variable heating power and emulates a borehole heat exchanger by a finite line-source model or a thermal resistance and capacity model. The possibly unknown parameters identified may comprise the thermal conductivity and volumetric heat capacity of the ground and grout, as well as the spacing between the pipes and the initial ground temperature. It is possible to integrate to the inversion the temperature measurements made at various depths in the fluid and grout and to take into account the fluid flow rate and the thermal capacity of the underground components.  

##Content of the repository
This repository contains the latest implementation of TRT-SInterp.  The folder entilted TRT-SInterp contains the source code, which is a collection of Matlab functions. The folder entilted Dataset contains an example file and a synthetic dataset.

##User manual
To learn how to use TRT-SInterp or learn about its theoretical foundations, please consult reference 1.

## Collaboration 
To suggest improvements, report a possible bug, or initiate collaboration, please, contact me at : philippe.pasquier@polymtl.ca or on [ResearchGate](https://www.researchgate.net/profile/Philippe_Pasquier2).

## Reference
Please, cite this work as: 

1. Pasquier, P., 2015. Stochastic interpretation of thermal response test with TRT-SInterp. Computers & Geosciences, 75, pp.73–87.


## Additional references
Additional information on TRT-SInterp can be found in the following references :
##
2. Pasquier, P. & Marcotte, D., 2013. Joint Use of Quasi-3D Response Model and Spectral Method to Simulate Borehole Heat Exchanger. Geothermics.
3. Pasquier, P. & Marcotte, D., 2012. Short-term simulation of ground heat exchanger with an improved TRCM. Renewable Energy, 46, pp.92–99.
4. Jacques, L., Pasquier, P., Marcotte, D.  Influence of Measurement and Model Error on Thermal Response Test Interpretation, in: Proceedings of the World Renewable Energy Congress 2014, London, United Kingdom. 2014.
5. Claesson, J. & Javed, S., 2011. An analytical method to calculate borehole fluid temperatures for timescales from minutes to decades. In ASHRAE Annual Conference. Montréal, Canada, p. 10.
6. Marcotte, D. & Pasquier, P. 2008.  On the estimation of thermal resistance in borehole conductivity test. Renewable Energy, vol. 33, p. 2407-2415. doi:10.1016/j.renene.2008.01.021
7. Marcotte, D. & Pasquier, P., 2008. Fast fluid and ground temperature computation for geothermal ground-loop heat exchanger systems. Geothermics, 37(6), pp.651–665.
8. Hellström, G., 1991. Ground Heat Storage. Thermal Analysis of Duct Storage Systems. Part I Theory. University of Lund,  Sweden.
9. Bennet, J., Claesson, J. & Hellström, G., 1987. Multipole Method to Compute  the Conductive Heat Flows to and between Pipes in a Composite Cylinder, Lund, Sweden: University of Lund, Department of Building Technology and Mathematical Physics.
