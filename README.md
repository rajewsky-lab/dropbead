# Dropseq - data analysis
This package contains classes and functions to analyze single cell 
sequencing data, like those generated through the Dropseq protocol.
The package is in essence tailored for the needs of Dropseq.

There are two main classes for the time being: 

* the `SingleSpeciesSample` contains objects with a `species1`,
`cells`, `genes` attributes and a `data.frame` representing DGE.
* the `MixedSpeciesSample` extends the `SingleSpeciesSample`
class and contains an additional `species2` attribute.

The following methods are contained in the package:
