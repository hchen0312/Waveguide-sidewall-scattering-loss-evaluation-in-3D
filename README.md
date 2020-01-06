# Waveguide-sidewall-scattering-loss-evaluation-in-3D
Computing the sidewall scattering loss of photonic waveguides using volume current method.

%%%%%% Important note %%%%%

Most importantly, to run this code on matlab, you should download the FDM solver from:

https://www.photonics.umd.edu/software/wgmodes/

Credit & Citation: Fallahkhair, Arman B., Kai S. Li, and Thomas E. Murphy. "Vector finite difference modesolver for anisotropic
dielectric waveguides." Journal of Lightwave Technology 26.11 (2008): 1423-1431.

We will call several functions from "wgmodes": wgmodes.m, postprocess.m, contourmode.m 

Alternatively, you can use your own FDM or FEM solver to compute the mode profile.

The obtained mode profile serves as the input for dyadic Green's function

%%%%%%

Some other useful references:

1. Chen, Hong, et al. "Study of crystalline defect induced optical scattering loss inside photonic waveguides in UV–visible spectral wavelengths using volume current method." Optics Express 27.12 (2019): 17262-17273.

2. Bauters, Jared F., et al. "Ultra-low-loss high-aspect-ratio Si 3 N 4 waveguides." Optics express 19.4 (2011): 3163-3174.

3. Barwicz, Tymon, and Hermann A. Haus. "Three-dimensional analysis of scattering losses due to sidewall roughness in microphotonic waveguides." Journal of Lightwave Technology 23.9 (2005): 2719.

The code is simplified for the purpose of better interpretation, alternatively, one can use it to compute scattering loss from threading dislocations as well, simply by removing some commented lines.

Running the code, you will obtain several plots:

Figure 1: Mode profile (Ex, Ey, Ez)

Figure 2: "Volume current" at the sidewall, whose radiation is resonsible for the scattering loss

Figure 3: Distributaion of Poynting vector in the far field in polar plot. Coordinate is defined in Ref. 1

neff indicates the effective index for the computed mode.

