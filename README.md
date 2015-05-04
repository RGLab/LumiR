LumiR: Luminex data analysis in R
=====

This package provide core functionality for reading and manipulating bead-level data, automated approaches for standard curve fitting 
and advanced bead-level analyses. 

# TODO

1. Does layout really have to be the same (e.g: not enough samples to fill last plate)
4. Add checks for concentration units in xml files
5. Allow one .lxd per plate. ->extract more info than just the analytes
6. Find other source for data (extdata is too big)

### Roxygenise
#### Methods
- show
- head
- pData
- fData
- exprs
- fit
- concentration
- subset
- formula
- getCoeffs
- merge
#### Functions
- results.curves.csv
- results.conc.csv

#### Done
- blum
- slum
- setup_templates
- slummarize
- read.experiment
- writeMBAA
- plot_layout
