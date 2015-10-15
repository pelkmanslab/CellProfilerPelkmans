# Tests

## To run tests
1. Set the environment variable *CELLPROFILER_PELKMANS_TESTS_DIR* to the location of test data e.g.
> S:\Data\Datasets\GitHubTests\CellProfilerPelkmans  (camelot-share-3)

2. Set the current matlab folder to *Scripts/*
3. Run the script *runAllTests.m*

## Types of Tests

### Handles-Injection 

CellProfiler modules accept a structure-array *handles* as an argument. Some processing is performed which alters *handles* and it is returned.

One approach is to *record* examples of this in-state and out-state. Then, in the future, run the module again on this in-state. The out-state should be identical to previously. 

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
Snapshots of execution state (e.g. serialized *handles* files) from matlab.
