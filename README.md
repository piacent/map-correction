# Map correction

Code for pixel-by-pixel application of map correction, to compute correction to `sc_integral` in the reco files. The new quantity is saved into a new leaf called `sc_integral_mapcorr`.

## To compile

```
g++ cygno-analyzer/Analyzer.cxx ApplyMapCorrection.cxx -o run.exe `root-config --libs --cflags` -lSpectrum
```

## Suggested Usage

```
./run.exe <path-to-input-recofile> <path-to-correction-map> <output-filename>
```

