CelltrackR README
================
Inge Wortel and Johannes Textor

<!-- README.md is generated from README.Rmd. Please edit that file -->

# CelltrackR

CelltrackR is designed to help with describing, visualizing, and
quantifying tracks of moving objects. Many common measures used in
physics and biology are implemented, such as mean square displacement
and autocorrelation. The package also provides a flexible function to
import tracks from text files.

The CelltrackR package has been developed as part of the MotilityLab
project. On the project website (motilitylab.net), a simple GUI frontend
to many functions in the package is implemented, and public datasets are
available to download and analyze. Indeed, the CelltrackR package was
formerly called “MotilityLab”, but we changed the name to make this
functionality easier to find.

See also: <https://ingewortel.github.io/celltrackR/>

# Installation

The latest development version of CelltrackR can be installed from
GitHub (note: this requires the `devtools` package):

``` r
   devtools::install_github( "ingewortel/celltrackR" )
```

Then the package can be loaded as usual:

``` r
    library(celltrackR)
```

# Examples

Tracks are organized as lists of matrices, and have S3 class `tracks`.
Three example datasets are provided with the package. A `plot` method is
implemented and can be used to visualize these datasets:

``` r
    plot( TCells, col=1 )
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

``` r
    plot( BCells, col=2 )
```

![](man/figures/README-unnamed-chunk-4-2.png)<!-- -->

To generate a simple mean square displacement plot for the example
dataset \`TCells’, use:

``` r
    msqd <- aggregate( TCells, squareDisplacement )
    plot( msqd, type='l' )
```

![](man/figures/README-unnamed-chunk-5-1.png)<!-- -->

This computes the squared displacement (MSD) for over all subtracks of
the dataset, and computes the average stratified by subtrack length. To
compute the MSD for non-overlapping subtracks only, use:

``` r
    msqd <- aggregate( TCells, squareDisplacement, max.overlap=0 )
    plot( msqd, type='l' )
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

MSD estimates can be biased in applications with a finite field of view
(such as microscopy), because slower objects remain in the field of view
for longer times. This can complicate comparisons between different
populations. We may thus wish to restrict our comparison to subtracks of
a certain (short) length. This can be done as follows:

``` r
    msqd.t <- aggregate( TCells, squareDisplacement, subtrack.length=1:5, max.overlap=0 )
    msqd.b <- aggregate( TCells, squareDisplacement, subtrack.length=1:5, max.overlap=0 )

    plot( msqd.t, type='l' )
    lines( msqd.b, col=2 )
```

![](man/figures/README-unnamed-chunk-7-1.png)<!-- -->

Another common measure to analyze tracks is the autocovariance function;
geometrically speaking, this is the dot product between pairs of pairs
of positions a fixed distance apart. For random walks, the
autocovariance decreases to 0 as the subtrack length increases. The
speed of convergence gives an indication of persistence of orientation,
a feature of many realistic objects. To compare persistence of the T
cells and B cells datasets, we can use:

``` r
    angle.t <- aggregate( TCells, overallDot )
    angle.b <- aggregate( BCells, overallDot )
    plot( angle.t, type='l' )
    lines( angle.b, col=2 )
```

![](man/figures/README-unnamed-chunk-8-1.png)<!-- -->

Many other ways to quantify tracks are implemented in the package and
described in the PDF documentation. For an overview of the available
commands, use the function

``` r
    help( package="celltrackR" )
```
