# Stacked-Films-Transfer-Matrix-Method




The scripts and functions are used to compute the following properties of a multi-layer thin-film stack. The fundamental function is "ObjectivePbI2_Au.m":
- electric field
- local absorption probability (dimensionless)
- Absorptivity of each layer ( [L]^(-1) dimension)
- Generation rate and spectral density ( [L]^(-3)[T]^(-1) and [L]^(-4)[T]^(-1))
- Reflection coefficient and power reflectivity at the front ( dimensionless)
- Transmission coeeficient and Power Transmissivity at the back (dimensionless)
- XYZ (or Yxy) chromaticity of the device in terms of reflected light

To run the function properly, it requires the following dataset:
- 'materials_nk.mat': the complex refractive indices of all involved layers
- 'AM15.mat': the AM1.5 GNI irradiance spectral density
- 'xyzMatrices.mat': the three CIE 1931 XYZ spectral color matching functions

The latter two are included. The first mat file needs to be extracted either from the previous datafrom or from the actual measurement of each layer material.

The calculation can run much faster if you enable the GPU computing in MATLAB if you have CUDA supported graphic cards. Simply replacing each matrix used in the computation by a gpuarray, e.g., "E=gpuarray(E);"; The final data can be transformed back for visulization by "E=gather(E);".

Make sure to created a seperate GPU-enabled version if you need it, so the CPU-only version can still be run on unsupported platforms.

-----------------------------------------------

To optimize for a given color, run the ParticleSwarmOptim.m. Be sure to replace the actual kernel function ('Objectivexxxx.m') by the version needed. The objective function is computed using ' ColorDelta.m'


