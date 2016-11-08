\name{motilityLab-package}
\alias{motilityLab-package}
\alias{motilityLab}
\docType{package}
\title{
MotilityLab: Quantitative Motion Analysis
}
\description{
Statistics to describe moving things
such as cells, bacteria or animals. 
Available measures include mean square displacement, 
confinement ratio, autocorrelation, asphericity,
turning angle, and fractal dimension.
}
\details{
\tabular{ll}{
Package: \tab motilityLab\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-04\cr
License: \tab GPL-2\cr
}
The MotilityLab package is designed for analyzing cell tracks acquired by
time-lapse microscopy (like those provided in the included datasets 
\code{\link{TCells}}, \code{\link{BCells}} and \code{\link{Neutrophils}}). 
But it can of course process any x-y-(z)-t data, and
we hope that it may be useful for other purposes as well.

The basic data structure that most functions in this package operate on is a set of 
\emph{tracks}. A track is a list of spatial coordinates that are recorded at 
\emph{fixed} time intervals; the function \code{\link{timeStep}} can be used to check
for fluctuations of the recording intervals. 
We expect tracks to be stored in a matrix (or data frame, but this is discouraged
for efficiency reasons) whose first column denotes a time interval 
(e.g. seconds elapsed since the beginning
of the experiment), and whose remaining columns denote a spatial coordinate. A set of
tracks is stored as a \code{\link{list}} with S3 class \code{tracks}. MotilityLab provides
some S3 methods for this class, which are explained in \code{\link{tracks}} as well as 
\code{\link{plot.tracks}}, \code{\link{sort.tracks}} and \code{\link{as.list.tracks}}.

A wide range of common track measures are included in the package. These are all functions
that take a single track as an input, and output one or several numbers.
For instance, the function \code{\link{speed}} estimates
the average instantaneous speed of the track by linear interpolation, and 
\code{\link{straightness}} computes the start-to-end distance divided by the trajectory
length (a number between 0 and 1, where 1 indicates a perfectly straight track).

MotilityLab is designed to support various flavors of track analysis that have been
suggested in the literature. The simplest kind is a \emph{track-based} analysis, where
we compute a single statistic for each track in a dataset (Beltman et al, 2009). Because
track sets are lists, this is achieved simply by using \code{lapply} or
\code{sapply} together with the track measure (see Examples).
 
In \emph{step-based} analyses (Beltman et al, 2009), we
chop each track up into segments of the same length and then apply our measures to those
segments. This can help to avoid biases that arise from variations
in track length (which are always present in cell tracking experiments).
In MotilityLab, step-based analyses are performed by using the \code{\link{subtracks}}
function. Often we want to perform such step-based analyses for all possible subtrack
lengths simultaneously, and plot the result as a function of the subtrack length;
a famous example is the \emph{mean square displacement plot}. This can be 
achieved by using the \code{\link{aggregate.tracks}} function, which has options
to control which subtrack lengths are considered and whether overlapping subtracks are 
considered.

In a \emph{staggered} \emph{staggered} analysis (Mokhtari et al, 2013), we analyse all
subtracks (of any length) of a single track, and typically plot the result as a matrix.
This can reveal dynamic patterns along a single track,
e.g. turning behaviour or local slowdowns. Staggered analyses can be performed using the
\code{\link{applyStaggered}} function.

Lastly, in addition to data analysis, the package contains some function to generate 
cell tracks by simulation. This is useful to develop and benchmark track analysis 
methodology (Textor et al, 2011), and for computational biology studies that try to 
extrapolate the long-term consequences of observed cell migration behaviour. Alongside
a simple uncorrelated random walk (\code{\link{brownianTrack}}), this package implements
a simulation model proposed by Beauchemin et al (2007) 
in the function \code{\link{beaucheminTrack}}. That model can also simulate 
directionally biased motion.

For a complete list of functions, use \code{library( help="MotilityLab" )}.
}
\author{
Katharina Dannenberg, Jeffrey Berry, Johannes Textor

Maintainer: Johannes Textor <johannes.textor@gmx.de>
}
\references{

Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009).
Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
789--798. doi:10.1038/nri2638

Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
and Marc Thilo Figge (2013), Automated Characterization and 
Parameter--Free Classification of Cell Tracks Based on Local Migration 
Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808

Johannes Textor, Antonio Peixoto, Sarah E. Henrickson, Mathieu
Sinn, Ulrich H. von Andrian and Juergen Westermann (2011),
Defining the Quantitative Limits of Intravital Two-Photon Lymphocyte Tracking.
\emph{PNAS} \bold{108}(30):12401--12406. doi:10.1073/pnas.1102288108

Catherine Beauchemin, Narendra M. Dixit and Alan S. Perelson (2007), Characterizing 
T cell movement within lymph nodes in the absence of antigen. \emph{Journal of Immunology}
\bold{178}(9), 5505-5512. doi:10.4049/jimmunol.178.9.5505
}
\keyword{cluster
spatial}
\examples{
## track-based speed comparison
boxplot(sapply( Neutrophils, straightness ), sapply( BCells, straightness ))

## step-based turning angle comparison
boxplot(sapply(subtracks(Neutrophils, 2), overallAngle), 
  sapply(subtracks(BCells, 2), overallAngle))

## mean square displacement plot; a step-based displacement analysis for all step lengths
plot(aggregate(TCells, squareDisplacement)[,"value"]) 

## 'staggered' analysis of displacement over whole track. Reveals that this track
## slows down near its beginning and near its end. 
filled.contour(applyStaggered(TCells[[4]], displacement, matrix=TRUE))

## a simple hierarchical clustering based on 2D asphericity

## tag track IDs so we can identify them later
names(TCells) <- paste0("T",names(TCells))
names(BCells) <- paste0("B",names(BCells))
names(Neutrophils) <- paste0("N",names(Neutrophils))

## project all tracks down to 2D
cells <-  projectDimensions(c(TCells,BCells,Neutrophils), c("x","y"))

## compute asphericity
asph <- lapply(cells, asphericity)

## plot clustering
plot(hclust(dist(asph)))
}
