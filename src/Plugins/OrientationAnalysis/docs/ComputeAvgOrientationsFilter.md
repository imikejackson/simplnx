# Compute Feature Average Orientations

## Group (Subgroup)

Statistics (Crystallography)

## Description

This **Filter** computes the average crystal orientation for each **Feature** (grain). Since each grain is made up of many **Cells** (voxels) that each have their own measured orientation, this filter combines those individual measurements into a single representative orientation per grain.

The average orientation is used by many downstream filters (e.g., misorientation calculations, Schmid factor, GBCD) and is one of the fundamental statistics computed during microstructure characterization.

Three averaging methods are available, and each can be independently enabled. Their results are stored in separate output arrays. The first method is the original DREAM3D Rodrigues average. The other two -- von Mises-Fisher and Watson -- are ports of the directional-statistics routines from the **EMsoft** package written by Dr. Marc De Graef's research group at Carnegie Mellon University. See *Origin of the vMF and Watson Methods* below.

### Method 1: Rodrigues Average (Original)

This is the original averaging algorithm. It determines the average orientation of each **Feature** by:

1. Gathering all **Elements** that belong to the **Feature**
2. Using the symmetry operators of the phase of the **Feature**, rotating the quaternion of the **Feature**'s first **Element** into the *Fundamental Zone* nearest to the origin
3. Rotating each subsequent **Element**'s quaternion (with same symmetry operators) looking for the quaternion closest to the current running average
4. Accumulating a running sum of the nearest quaternions
5. Dividing the accumulated quaternion sum by the count and normalizing to produce the average

The process of finding the nearest quaternion in Step 3 accounts for the periodicity of orientation space, which would cause problems in the averaging if all quaternions were forced to be rotated into the same *Fundamental Zone*. The quaternions can be averaged with a simple summation-based average because quaternion space is not distorted like Euler space.

**Outputs:** Average Quaternions, Average Euler Angles (Bunge convention Z-X-Z)

### Method 2: Von Mises-Fisher (vMF) Average

The von Mises-Fisher distribution is a probability distribution on the surface of a unit hypersphere in *p*-dimensional space. For orientation averaging, the relevant case is the unit quaternion sphere (*p* = 4). The vMF distribution is parameterized by:

- **mu** (mean direction): A unit quaternion representing the central tendency of the distribution. This is the estimated average orientation.
- **kappa** (concentration parameter): A non-negative scalar that characterizes how tightly the orientations are clustered around the mean. Intuitively, kappa is to the vMF distribution what the full-width-at-half-maximum (FWHM) is to a Gaussian distribution: it is a measure of how narrow or tight the distribution is. Higher kappa values indicate tighter clustering (less spread); kappa = 0 corresponds to a uniform distribution on the sphere.

The vMF probability density for a unit vector **x** given mean direction **mu** and concentration **kappa** is proportional to exp(kappa * mu^T * x). This makes it the spherical analogue of the Gaussian distribution on a flat space.

A crystal orientation is only defined up to the symmetry operators of its Laue class, so the filter does not fit a plain vMF distribution. It fits the **modified (symmetry-group-invariant) von Mises-Fisher distribution**, in which the density is an equal-weight mixture of one vMF component per proper rotation operator of the phase's Laue class (24 components for cubic, 1 for triclinic). This is the formulation derived by Chen and co-workers [2] and implemented in EMsoft.

The filter estimates the vMF parameters using an **Expectation-Maximization (EM)** algorithm. All element quaternions belonging to a feature are first reduced to the *Fundamental Zone* using the crystal symmetry operators. The EM procedure then iteratively refines the estimates of **mu** (the average orientation quaternion) and **kappa** (the concentration).

**Outputs:** Average Quaternions, Average Euler Angles (Bunge convention Z-X-Z), Kappa Values

### Method 3: Watson Average

The Watson distribution is a probability distribution on the unit sphere that is **antipodally symmetric**, meaning it treats **x** and **-x** as equivalent. This property makes it particularly well-suited for orientation data represented as quaternions, since quaternions **q** and **-q** represent the same physical rotation.

The Watson distribution is parameterized by:

- **mu** (mean axis): A unit quaternion representing the principal axis of the distribution. This is the estimated average orientation.
- **kappa** (concentration parameter): A scalar that controls the concentration of the distribution around the mean axis. As with the vMF distribution, kappa is analogous to the full-width-at-half-maximum (FWHM) of a Gaussian: it measures how narrow or tight the distribution is. For positive kappa, the distribution is bipolar (concentrated around +/-mu); for negative kappa, it is girdle-shaped (concentrated in the great circle perpendicular to mu).

The Watson probability density for a unit vector **x** is proportional to exp(kappa * (mu^T * x)^2). The key difference from the von Mises-Fisher distribution is the squared dot product, which enforces the antipodal symmetry.

As with the vMF method, the filter fits the **modified (symmetry-group-invariant) axial Watson distribution** -- an equal-weight mixture of one Watson component per proper rotation operator of the phase's Laue class -- rather than a plain Watson distribution [2].

Like the vMF method, the filter estimates Watson parameters using an **Expectation-Maximization (EM)** algorithm operating on fundamental-zone-reduced quaternions.

**Outputs:** Average Quaternions, Average Euler Angles (Bunge convention Z-X-Z), Kappa Values

### Origin of the vMF and Watson Methods

The von Mises-Fisher and Watson averaging methods are ports of the **directional statistics** routines from the [EMsoftOO package] (https://github.com/EMsoft-org/EMsoftOO/blob/develop/Source/EMsoftOOLib/mod_DIC.f90), written by **Dr. Marc De Graef's research group at Carnegie Mellon University**, with the original Expectation-Maximization implementation contributed by Yu-Hui Chen (University of Michigan). The statistical formulation is published in Chen *et al.* [1] [2].

EMsoft is distributed under a BSD 3-Clause license, Copyright (c) 2014-2022 Marc De Graef Research Group / Carnegie Mellon University.

### Hard-Coded Algorithm Parameters

The following parameters are currently hard-coded in the implementation and are not user-configurable:

| Parameter | EMsoft name | Value | Description |
|-----------|-------------|-------|-------------|
| **Random Seed** | `seed` | 43514 | Seed for the pseudo-random number generator that draws the starting guess for **mu** at the beginning of each EM restart. Because the seed is fixed, the vMF and Watson results are deterministic and reproducible from run to run. |
| **EM Restarts** | `Num_of_init` | 5 (dimensionless count) | Number of independent EM runs, each begun from a different random unit quaternion. Expectation-Maximization can converge to a local maximum of the likelihood, so the algorithm restarts several times and keeps the solution with the highest likelihood. |
| **EM Iterations** | `Num_of_iterations` | 10 (dimensionless count) | Maximum number of Expectation-Maximization iterations performed within each restart. The loop exits early once the Q-function changes by less than 0.01 between successive iterations. |

These values may be exposed as user-configurable parameters in a future release.

### Special Cases

- **Features with a single element:** For the vMF and Watson methods, if a feature contains only one element orientation, the EM algorithm is skipped entirely and the single quaternion is used directly as the average. The kappa value is set to 0 in this case.
- **Features with zero elements:** Features with no elements (phase <= 0 for all voxels) will have their output arrays initialized to NaN (for vMF/Watson) or identity quaternion / zero Euler angles (for Rodrigues).
- **Phase indexing:** The filter requires that phase values be > 0 for elements to be included in the averaging. Phase index 0 is reserved for "Unknown" in the Crystal Structures array and is always skipped. This applies identically to all three methods.
- **Invalid phases and crystal structures:** Elements whose phase value lies outside the range of the Crystal Structures array, and elements or features whose crystal structure value is not a supported Laue class (for example 999 = Unknown), are **excluded** from the averaging. The filter emits a warning (-54672 for out-of-range phases, -54671 for unknown crystal structures) reporting how many were dropped — the drop is never silent. A feature whose elements are all excluded finalizes to the identity quaternion (Rodrigues) or NaN (vMF/Watson).
- **Multi-phase features:** The vMF/Watson methods use a single crystal structure per feature, taken from the phase of the feature's highest-index element; the Rodrigues method uses each element's own phase. Features are normally single-phase, so this distinction rarely matters.
- **No method enabled:** If none of the three averaging methods is enabled the filter fails in preflight with error -54673.

### Required Input Sources

- **Cell Quaternions** -- typically read from EBSD data via [Read H5EBSD](ReadH5EbsdFilter.md), [Read CTF Data](ReadCtfDataFilter.md), or [Read ANG Data](ReadAngDataFilter.md); can also be produced from Euler angles by [Convert Orientations](ConvertOrientationsFilter.md).
- **Cell Feature Ids** -- produced by a segmentation filter such as [Segment Features (Misorientation)](EBSDSegmentFeaturesFilter.md).
- **Cell Phases** -- typically read from EBSD data alongside the quaternions.
- **Crystal Structures** -- ensemble-level array read from EBSD data or created by [Create Ensemble Info](CreateEnsembleInfoFilter.md).

% Auto generated parameter table will be inserted here

## References

[1] Y.H. Chen, S.U. Park, D. Wei, G. Newstadt, M.A. Jackson, J.P. Simmons, M. De Graef, and A.O. Hero. A Dictionary Approach to Electron Backscatter Diffraction Indexing. *Microscopy and Microanalysis*, 21(3), 739-752 (2015). DOI: [10.1017/S1431927615000756](https://doi.org/10.1017/S1431927615000756)

[2] Y.H. Chen, D. Wei, G. Newstadt, M. De Graef, J.P. Simmons, and A.O. Hero. Parameter Estimation in Spherical Symmetry Groups. *IEEE Signal Processing Letters*, 22(8), 1152-1155 (2015). DOI: [10.1109/LSP.2014.2387206](https://doi.org/10.1109/LSP.2014.2387206)

[3] EMsoft, Marc De Graef Research Group, Carnegie Mellon University. Module `dictmod` (`Source/EMsoftLib/dictmod.f90`). [https://github.com/EMsoft-org/EMsoft](https://github.com/EMsoft-org/EMsoft)

[4] K.V. Mardia and P.E. Jupp. *Directional Statistics*. Wiley, Chichester (2000). Background reference for the von Mises-Fisher and Watson distributions.

## Example Pipelines

+ (02) Small IN100 Full Reconstruction
+ (05) SmallIN100 Crystallographic Statistics
+ (06) SmallIN100 Synthetic

## License & Copyright

Please see the description file distributed with this **Plugin**

## DREAM3D-NX Help

If you need help, need to file a bug report or want to request a new feature, please head over to the [DREAM3DNX-Issues](https://github.com/BlueQuartzSoftware/DREAM3DNX-Issues/discussions) GitHub site where the community of DREAM3D-NX users can help answer your questions.
