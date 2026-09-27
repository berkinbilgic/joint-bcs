# Multi-Contrast Reconstruction with Bayesian Compressed Sensing

MATLAB code and numerical phantoms for joint Bayesian compressed sensing (BCS) reconstruction of multi-contrast images from undersampled k-space.

- `script_Bayesian_CS_realValued.m`: joint BCS reconstruction assuming real-valued images
- `script_Bayesian_CS_complexValued.m`: complex-valued joint and non-joint reconstructions

## Usage

Run either script from the repository folder in MATLAB.

## Third-party code

`utils/` includes functions from Michael Lustig's SparseMRI code and Gabriel Peyre's toolboxes (copyright notices retained in the files). The CurveLab curvelet transform files from the original distribution (`fdct_wrapping*`, `ifdct_wrapping`, `fdct2vec`, `vec2fdct`) are not included, because CurveLab has its own license; neither script uses them.

## Reference

B Bilgic, VK Goyal, E Adalsteinsson. Multi-contrast Reconstruction With Bayesian Compressed Sensing. *Magnetic Resonance in Medicine*, 2011. [[PDF]](https://www.martinos.org/~berkin/Bilgic_2011_BCS_MRM.pdf)

Contact: berkin AT nmr.mgh.harvard.edu
