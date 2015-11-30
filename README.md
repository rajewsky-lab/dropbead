# Dropseq - data analysis
This package contains classes and functions to analyze single cell 
sequencing data, like those generated through the Dropseq protocol.
The package is in essence tailored for the needs of Dropseq.

There are three main classes: 

* the `DigitalGeneExpressionMatrix` class is basically the DGE
and contains operations on DGE
* the `SingleSpeciesSample` contains objects with a `species1`
attribute and a `DigitalGeneExpressionMatrix`
* the `MixedSpeciesSample` extends the `SingleSpeciesSample`
class and contains an additional `species2` attribute.

The following methods are contained in the package:
