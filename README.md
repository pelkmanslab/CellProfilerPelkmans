CellProfilerPelkmans
====================

A customized-version of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
* many additional custom modules

## How to install?

CellProfilerPelkmans is part of each [iBRAIN_UZH release](https://github.com/pelkmanslab/iBRAIN_UZH/releases) in *iBRAIN\tools\CellProfilerPelkmans*.

Please follow the instructions in the [iBRAIN_UZH User Guide](https://github.com/pelkmanslab/iBRAIN_UZH/blob/master/doc/USER_GUIDE.md)

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

Reference documentation is found in several locations:
* [Modules API Documentation](http://jenkins.pelkmanslab.org/job/CellProfilerPelkmans_Master/CellProfilerPelkmans_API_Documentation/workspace/Modules/index.html). It is updated automatically daily, and contails all modules.
* [Core Module help](https://github.com/pelkmanslab/iBRAIN_BRUTUS/wiki/iBRAIN_BRUTUS-core-module-help). A subset of *core* modules
* [CellProfiler manual](http://cellprofiler.org/linked_files/Documentation/cp1_manual_9717.pdf). Standard modules only.
* [CellProfilerPelkmans Wiki](https://github.com/pelkmanslab/CellProfilerPelkmans/wiki). Custom modules only.
* [List of Standard Modules](doc/LIST_OF_STANDARD_MODULES.md). Standard modules only.

Further documentation can be found in the actual Matlab functions (.m files). This information can also be queried from the GUI using the `?` button.
