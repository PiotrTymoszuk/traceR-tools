# traceR-tools
Development of QC wrappers for cell - cell pair statistics and doublet identification

Tools for calculation of cell - cell position correlations (Spearman) and distances wrapped around the object 'traces' provided by 'tracer' package (https://cran.r-project.org/web/packages/tracer/index.html). A real experiment example provided. 

Cell doublet identification based on the assumption: they should stay near each other (low mean distance and low SD of the cell - cell distance).

Pending: identification of temporal cell association and track crossing based on cell - cell distance minima.

Will turn into a package in the future:)
