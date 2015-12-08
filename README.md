# Dropseq - data analysis
This package contains classes and functions to analyze single cell 
sequencing data, like those generated through the Dropseq protocol.
The package is in essence tailored for the needs of Dropseq.

There are two main classes: 
* the `SingleSpeciesSample` contains objects with a `species1`,
`cells`, `genes` attributes and a `data.frame` representing DGE.
* the `MixedSpeciesSample` extends the `SingleSpeciesSample`
class and contains an additional `species2` attribute.

Methods contained in the `SingleSpeciesSample` class:
* `computeGenesPerCell`.
* `computeTranscriptsPerCell`.
* `geneExpressionMean`.
* `geneExpressionDispersion`.
* `geneExpressionVariability`.
* `removeLowQualityGenes`, removes genes which are not expressed in
a minimum number of cells.
* `removeLowQualityCells`, removes cells which do not express a minimum
number of genes.
* `listCellsToCollapse`, returns a list with pairs of cells sharing
11 bases in their barcodes.
* `collapseCellsByBarcodes`.

Methods contained in the `MixedSpeciesSample` class (in addition
to the ones that exist through polymorphism):
* `splitDgeByGenesOfSpecies`, returns a list of two `data.frames`
corresponding to genes of each species. All cells are included.
* `classifyCellsAndDoublets`, returns a `data.frame` with cells,
counts of transcripts per species and characterization of cell.
* `splitDgeByGenesAndCellsOfSpecies`, returns a list of two 
`data.frames` corresponding to genes and cells of each species.
* `splitMixedSpeciesSampleToSingleSpecies`, returns a list of 
two `SingleSpeciesSample` objects separating the species.
