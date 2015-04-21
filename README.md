CellProfilerPelkmans
====================

A fork of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
* refreshed GUI 
* many additional custom modules


## Installation ##

First, you need to clone this repository:

```{bash}
git clone git@github.com:pelkmanslab/CellProfilerPelkmans.git ~/cellprofiler
```

In order to be able to start CellProfiler from Matlab, you have to put the location of the local copy of this repository on the Matlab path.
To this end, add the following line to your `startup.m` file:

```{matlab}
addpath(genpath('~/cellprofiler'))
```

Note that this file has to reside in your *initial working folder*. For more information see [Matlab startup](http://ch.mathworks.com/help/matlab/ref/startup.html).

## Dependencies ##

Some CP modules depend on Matlab code that resides outside of this repository. This code is located at [PelkmansLibrary](https://github.com/pelkmanslab/PelkmansLibrary).

To use *CellProfilerPelkmans* you thus also have to put the location of the local copy of the *iBRAINdependencies* repository on the Matlab path.

To this end, first clone the repository:

```{bash}
git clone git@github.com:pelkmanslab/iBRAINDependencies.git ~/pelkmanslibrary
```

Then add the following line to your `startup.m` file:

```{matlab}
addpath(genpath('~/pelkmanslibrary/matlab'))
```

On unix systems you can alternatively define the `MATLABPATH` environment variable. To this end, add the following line to your `.bash_profile` file:

```{bash}
export MATLABPATH=$MATLABPATH:$HOME/pelkmanslibrary/matlab
```


## How to ##

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
