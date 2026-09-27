# Multi-Contrast Reconstruction with Bayesian Compressed Sensing

MATLAB code and numerical phantoms for joint Bayesian compressed sensing (BCS) reconstruction of multi-contrast images from undersampled k-space.

- `script_Bayesian_CS_realValued.m`: joint BCS reconstruction assuming real-valued images
- `script_Bayesian_CS_complexValued.m`: complex-valued joint and non-joint reconstructions

## Usage

Run either script from the repository folder in MATLAB.

## Known issues

**Linux and macOS:** `script_Bayesian_CS_realValued.m` adds its helper folder with a Windows-style path (`addpath utils\`), which does not resolve, so the script stops with `Unrecognized function or variable 'tile'`. Add the folder first:

```matlab
addpath utils
script_Bayesian_CS_realValued
```

**`script_Bayesian_CS_complexValued.m`:** the separate (non-joint) reconstruction on lines 91-97 assumes a scalar image size, while the script sets `im_size` to a 2-element vector, so it stops at line 91 with `Incorrect dimensions for raising a matrix to a power`. Changing these four lines fixes it:

- lines 91-92: call `mt_CSfft2_rect` instead of `mt_CSfft2` (as the joint reconstruction on lines 107-108 already does)
- lines 96-97: reshape with `[im_size,2*L]` instead of `[im_size,im_size,2*L]`

With these changes both scripts run to completion (tested with MATLAB R2026a).

## Third-party code

`utils/` includes functions from Michael Lustig's SparseMRI code and Gabriel Peyre's toolboxes (copyright notices retained in the files). The CurveLab curvelet transform files from the original distribution (`fdct_wrapping*`, `ifdct_wrapping`, `fdct2vec`, `vec2fdct`) are not included, because CurveLab has its own license; neither script uses them.

## Reference

B Bilgic, VK Goyal, E Adalsteinsson. Multi-contrast Reconstruction With Bayesian Compressed Sensing. *Magnetic Resonance in Medicine*, 2011. [[PDF]](https://www.martinos.org/~berkin/Bilgic_2011_BCS_MRM.pdf)

Contact: berkin AT nmr.mgh.harvard.edu
