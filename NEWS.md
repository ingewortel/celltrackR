# celltrackR 1.0.0

Potentially breaking changes
- replaced the in-package datasets "TCells", "BCells", "Neutrophils" with new data;
this may break some examples using the old datasets because the new datasets have tracks
with different IDs, and because they are now 2D projections rather than 3D tracks.
Details on the data will be available in the celltrackR publication; raw data are 
available from inst/extdata (see the vignettes on QC for how to access those).
- changed the behavior of some of the angle functions like angleToPoint, angleToDir,
and angleToPlane. In the old version, if the supplied track was more than one step long,
only the angle from the first step to the reference was considered. This was not very 
compatible with how the other track measures work. The new version therefore computes
angles with the overall displacement vector of the supplied track. This should not 
alter conclusions too much, but it will change the numeric outputs of some of these
functions when applied to longer tracks.

Documentation changes:
- All the vignettes were updated to work with the new datasets
- The vignette on QC contains a new section on filtering out non-motile tracks
- A new vignette explicitly describes the preprocessing choices made for the package
datasets
- The vignette on track simulation now contains examples of how to fit models based
on the mean squared displacement curve.

Minor changes:
- Moved pracma from 'suggests' to 'imports' as it is now used by multiple functions
that work with angles
- as.tracks.data.frame now automatically handles a dataframe with 2D tracks if the
pos.columns argument is left unspecified; earlier, this would have resulted in an error.
