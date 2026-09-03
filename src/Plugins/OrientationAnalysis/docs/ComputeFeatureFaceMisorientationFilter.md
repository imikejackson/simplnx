# Compute Feature Face Misorientation (Face)

## Group (Subgroup)

Processing (Crystallography)

## Description

This **Filter** computes, for every **Triangle** in a **Triangle Geometry** surface mesh, the misorientation angle between the two **Features** that lie on either side of that **Triangle**. Misorientation is the angular difference between the crystal orientations of two grains -- it describes how much one grain's crystal lattice is rotated relative to another's.

![Fig. 1: The misorientation between two grains is the rotation relating their crystal orientations; its angle is the misorientation angle. Because a crystal has symmetry-equivalent orientations, the smallest equivalent angle — the disorientation — is what the filter reports.](Images/Misorientation_Concept.png)

The output is a **single scalar value per Triangle**: the misorientation angle in **degrees**, stored as a 1-component `float32` array. The array is created in the same **Face Data** attribute matrix as the selected **Face Labels** array, and is named `Face Misorientations` by default.

The angle is symmetry-reduced. Because a crystal has symmetry-equivalent orientations, the filter reports the smallest equivalent angle -- the *disorientation* -- using the Laue class of the phase the two **Features** belong to. This filter works for all valid crystal structures / Laue classes.

The axis of the misorientation is **not** computed or stored -- only the angle.

### Note

A value of NaN is stored for a **Triangle** in any of the following cases:

- Either **Face Label** is not a valid **Feature** id (that is, it is `0` or negative). These are the triangles on the outer surface of the volume, which have a **Feature** on one side only.
- The two **Features** belong to **different phases**. Misorientation between different crystal systems is not physically meaningful.
- The phase's **Crystal Structure** value does not correspond to a known Laue class.

### Required Input Sources

- **Face Labels** -- a 2-component `int32` face array naming the **Feature** on either side of each **Triangle**, produced by a surface meshing filter such as [Quick Surface Mesh](../SimplnxCore/QuickSurfaceMeshFilter.md).
- **Average Quaternions** -- produced by [Compute Average Orientations](ComputeAvgOrientationsFilter.md).
- **Feature Phases** -- produced by [Compute Feature Phases](../SimplnxCore/ComputeFeaturePhasesFilter.md).
- **Crystal Structures** -- ensemble-level array read from EBSD data or created by [Create Ensemble Info](CreateEnsembleInfoFilter.md).

### Change in Output Format

Earlier versions of this filter wrote a 3-component array of "axis \* angle" values. It now writes the single misorientation angle described above. A pipeline saved against the older version must have this filter's output parameter updated by hand; no automatic conversion is performed.

% Auto generated parameter table will be inserted here

## Example Pipelines

+ (07) Small IN100 Mesh Statistics

## License & Copyright

Please see the description file distributed with this **Plugin**

## DREAM3D-NX Help

If you need help, need to file a bug report or want to request a new feature, please head over to the [DREAM3DNX-Issues](https://github.com/BlueQuartzSoftware/DREAM3DNX-Issues/discussions) GitHub site where the community of DREAM3D-NX users can help answer your questions.
