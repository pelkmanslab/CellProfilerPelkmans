# Tests

## To run tests
1. Set the environment variable *CELLPROFILER_PELKMANS_TESTS_DIR* to the location of test data e.g.
> S:\Data\Datasets\GitHubTests\CellProfilerPelkmans

2. Set the current matlab folder to *Scripts/*
3. Run the script *runAllTests.m*

## Test data

The test-data folder is structured as follows:

### *Input/*
Raw-inputs for pipelines (images, illumination correction files etc.). Sub-directories *A, B, C, etc.* each represent different sets of raw-inputs.

### *State/*
Snapshots of execution state (e.g. serialized *handles* files) from matlab. They are grouped:

1. The purpose of the test (e.g. a module name)
2. Then by *A, B, C, etc.* referring to the associated input data
3. Then by a number in the filename *1_handles_in.mat*, as there might be several cases. Please preserve the filename pattern.

Each of these folder has a file *pathOld.txt* that records the images-directory used when the snapshops were made. This will be recorded within the *handles* data structures. It is used for search-and-replace operations during the unit tests, to update the paths to a new location.
