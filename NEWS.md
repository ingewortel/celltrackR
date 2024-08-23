# celltrackR 1.2.1

Improved handling of online resources in tests and examples for immunemap integration
functions in accordance with CRAN policy.

### Potentially breaking changes: 

None.

### Minor changes:
- (immunemap-integration.R => parse.immap.json ): now uses jsonlite instead of rjson 
	because rjson requires R 4.4.1.
- (immunemap-integration.R => parse.immap.json ): improved handling of cases where the 
	internet or immunemap website is down.
- (immunemap-integration.R ): dontrun on examples to avoid errors on CRAN servers due to
	internet problems.
- (testImmuneMapIntegration.R) : added a skip_if_offline for same reasons.
- DESCRIPTION: jsonlite instead of rjson in suggests


### Documentation changes:

Minor changes:
- Updated outdated @docType package to new roxygen recommendation


# celltrackR 1.2.0

Integration with immunemap.org. 

### Potentially breaking changes: 

None.

### Minor changes:
New functionalities for integration with immunemap.org:
- (immunemap-integration.R => read.immap.json ): new function, reads data from an 
	immunemap json file (or an url pointing to such a file). This returns a list with
	both track object(s) and metadata; see documentation ?read.immap.json. 
- (immunemap-integration.R) : see above; file contains also several helper functions.
- tests added for these reading functions as well.
- (inst/extdata/immunemap.json) : new example file showing the format for immunemap json.


### Documentation changes:

Minor changes:
- including documentation for new functionalities for integration with immunemap.org;
- Updated package documentation after breaking change in Roxygen upon CRAN request.
- Updated the documentation for the included datasets, which had not yet been updated 
with the new versions of the data included since package v 1.0.0.


# celltrackR 1.1.0

Significant speed-up in some of the angle analysis functions.

### Potentially breaking changes: 

None.


### Minor changes:
Most of the changes are internal only to speed some of the computations (especially
computing the pairwise angles and distances in analyzeCellPairs and analyzeStepPairs).
To accomplish this, some functions now have additional input arguments or can handle
inputs in an additional format (old input formats are also still okay):
- (angle-functions.R => pairsByTime ) 
	new function. At each time point in the dataset, looks for pairs of cells that co-occur
	at that time and computes their distance to each other. Outputs pairs, times of co-occurrence,
	and distances in a dataframe. Internally, this function uses the "dist" function to 
	compute pairwise distances rapidly; several of the other distance/angle functions 
	now use pairsByTime() internally for efficiency reasons.
- (angle-functions.R => angleSteps ) 
	Instead of only computing a single angle for two trackIDs at one timePoint,
	this function can now also compute multiple angles at the same time for efficiency reasons.
	Instead of a vector with two trackIDs, it can now also take a two-column matrix with
	one pair of IDs per row; the "t" argument must then also be a numeric vector instead
	of a single value. See the documentation under ?angleSteps.
- (angle-functions.R => angleCells )
	Input now allows to compute multiple distances at once; see angleSteps changes above.
- (angle-functions.R => distanceSteps )
	Input now allows to compute multiple distances at once; see angleSteps changes above.
- (angle-functions.R => distanceCells )
	Input now allows to compute multiple distances at once; see angleSteps changes above.
	Also added a 'quietly' argument to suppress warnings for NA distances.
- (angle-functions.R => analyzeStepPairs )
	Now has an argument "searchRadius" to only return pairs that start within a given
	distance from each other. 
- (angle-functions.R => analyzeCellPairs )
	Now has an argument "searchRadius" to only return pairs that start within a given
	distance from each other. 


### Bugfixes:
- (functions.R => timePoints )
	When supplied with a tracks objects where tracks had equal length, the 
	function erroneously returned a matrix rather than a single vector of timepoints
	(although this would have almost never happened). This is now fixed.
	Also ensures the list of timePoints is sorted.

### Documentation changes:

None, other than updating the documentation of the functions that can now take multiple
inputs simultaneously or have additional input arguments (see Minor changes).



# celltrackR 1.0.0

### Potentially breaking changes
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

### Documentation changes:
- All the vignettes were updated to work with the new datasets
- The vignette on QC contains a new section on filtering out non-motile tracks
- A new vignette explicitly describes the preprocessing choices made for the package
datasets
- The vignette on track simulation now contains examples of how to fit models based
on the mean squared displacement curve.

### Minor changes:
- Moved pracma from 'suggests' to 'imports' as it is now used by multiple functions
that work with angles
- as.tracks.data.frame now automatically handles a dataframe with 2D tracks if the
pos.columns argument is left unspecified; earlier, this would have resulted in an error.
