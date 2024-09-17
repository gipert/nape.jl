# Not a Precision Experiment (NAPE)

Global bayesian analysis of Ge-76 neutrinoless double beta decay data.
Currently supports data sets from the MAJORANA DEMONSTRATOR, GERDA and
LEGEND-200.

Quick start:
```
$ julia main.jl
```

## Useful pointers

- https://legend-exp.org
- https://bat.github.io/BAT.jl
- https://www.mpi-hd.mpg.de/gerda
- https://www.npl.washington.edu/majorana/majorana-experiment

## Notes

- fit window: [1930, 2190] keV, excluding ± 5 keV around 2104 keV from Tl-208
  and 2119 keV from Bi-212 -> 240 keV wide window
