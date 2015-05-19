# Tests

## To run tests
1. Set the environment variable *CELLPROFILER_PELKMANS_TESTS_DIR* to the location of test data e.g.
> S:\Data\Datasets\GitHubTests\CellProfilerPelkmans

2. Set the current matlab folder to *Scripts/*
3. Run the script *runAllTests.m*

## Types of Tests

### Handles-Injection 

CellProfiler modules accept a cell-array *handles* as an argument. Some processing is performed which alters *handles* and it is returned.

One approach is to *record* examples of this in-state and out-state. Then, in the future, run the module again on this in-state. The out-state should be identical to previously. 

In practice, some rewriting of paths needs to occur, as filesystem paths will differ in the *handles* cell-array between environments:

* **recording environment** - where in-state and out-state are recorded (one particular personal computer)
* **testing environment** - for testing in the future (on Jenkins server, varying personal computers) 

The value of the test:

* Provides a sanity-check that each [tested] module continuously works. But only on particular data and arguments.

Please see [Creating a Handles-Injection Test](CREATING_HANDLES_INJECTION_TEST.md).

## Test data

The test-data folder is structured as follows:

### *Input/*
Raw-inputs for pipelines (images, illumination correction files etc.). Sub-directories *A, B, C, etc.* each represent different sets of raw-inputs.

Within each folder, follow a typical pipeline structure:
* *TIFF* for image folders
* *BATCH* for measurement inputs and outputs
* Store pipelines in root folder.

### *State/*
Snapshots of execution state (e.g. serialized *handles* files) from matlab. They are grouped:

1. The purpose of the test (e.g. a module name)
2. Then by *A, B, C, etc.* referring to the associated input data
3. Then by a number in the filename *1\_handles\_in.mat*, as there might be several cases. Please preserve the filename pattern.

Each of these folder has a file *pathOld.txt* that records the images-directory used when the snapshops were made. This will be recorded within the *handles* data structures. It is used for search-and-replace operations during the unit tests, to update the paths to a new location.