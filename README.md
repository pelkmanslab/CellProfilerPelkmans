CellProfilerPelkmans
====================

A fork of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
* many additional custom modules

## How to install?

CellProfilerPelkmans is part of each [iBRAIN_UZH release](https://github.com/pelkmanslab/iBRAIN_UZH/releases). Please follow the instructions in the [iBRAIN_UZH User Guide](https://github.com/pelkmanslab/iBRAIN_UZH/blob/master/USER_GUIDE.md)

## Dependencies ##

Some CellProfilerPelkmans modules depend on Matlab code that resides outside of this repository - specifically in [PelkmansLibrary](https://github.com/pelkmanslab/PelkmansLibrary). This is also part of each iBRAIN_UZH release.

## After installation

Now you are ready to go.

Open your Matlab application and type:
```{matlab}
CellProfiler
```

This will start the CellProfiler program and open the GUI window. 


## Modules ##

Below you find a list of available standard Cell Profiler (CP) modules.

For documentation on these standard modules see [CellProfiler manual](http://cellprofiler.org/linked_files/Documentation/cp1_manual_9717.pdf).

For documentation on custom modules please refer to [CellProfilerPelkmans wiki](https://github.com/pelkmanslab/CellProfilerPelkmans/wiki).

Further documentation can be found in the actual Matlab functions (.m files). This information can also be queried from the GUI using the `?` button.

### Standard modules ###

##### File Processing #####

* ExportToExcel
* LoadSingleImage
* SaveImages (+)

##### Image Processing #####

* ApplyThreshold
* Combine
* Crop
* RescaleIntensities
* Resize
* Smooth
* Subtract
* SubtractBackground

##### Object Processing #####

* ExandOrShrink
* FilterByObjectMeasurement
* IdentifyPrimAutomatic
* IdentifyPrimLoG
* IdentifyPrimManual
* IdentifySecondary
* IdentifyTertiarySubregion
* Relate
* MergeAndRelateChildren

##### Measurements #####

* MeasureCorrelation
* MeasureImageGranularity
* MeasureImageIntensity
* MeasureObjectAreaShape
* MeasureObjectIntensity
* MeasureObjectNeighbors
* MeasureRadialDistribution
* MeasureTexture
* MeasureGPperSingleCell
* MeasureObjectColocalisation


##### Other #####

* SpeedUpCellprofiler


(+) indicates standard modules with custom improvements (e.g.: additional optional input arguments).
