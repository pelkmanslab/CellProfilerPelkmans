# Creating a Handles-Injection Test

A handles-injection test is a particular type of check performed on CellProfiler modules. 

CellProfiler modules accept a structure-array *handles* as an argument. Some processing is performed which alters *handles* and it is returned.
 
This type of test verifies that the behaviour of a module remains consistent by recording and comparing examples of *handles*.

Basically, the check performed is:

> module\_function( *recorded\_handles\_in\_state* ) = *recorded\_handles\_out\_state*

There are three steps necessary:

1. **Extract test data**. Find an example pipeline which is known to function correctly (for the module of interest). Select a small subset of images.
2. **Record handles-state**. Run this pipeline locally on your PC. Record *handles* when entering and leaving the module.
3. **Write a verification test**. Verify that the module produces the same output in future.

### Rewriting paths
In practice, the some rewriting of paths needs to occur, as the filesystem paths stored in *handles* will differ between environments:

* **recording environment** - where in-state and out-state are recorded (one particular personal computer)
* **testing environment** - for testing in the future (on Jenkins server, varying personal computers) 

So the **recording-path** is saved with the in-state and out-state. Paths are then rewritten as needed with a *search-and-replace* operation.

## Step 1: Extract test data 

1. Find an example pipeline that functions correctly. Make a new directory in:

	> **camelot-share-3**/Data/Datasets/GitHubTests/CellProfilerPelkmans/Input

	Follow the directory naming of A, B, C, etc. Pick the next available letter.

2. Copy the pipeline into this directory.

3. Make a subfolder *TIFF/*. Copy some images there. Ensure all channels exist.
 
4. Make a subfolder *BATCH/*. If necessary, copy illumination correction files (e.g. *Measurements\_batch\_illcor\_channel002\_zstack000.mat*) there.  

This test-data can be reused for testing multiple modules.

## Step 2: Record handles-state

1. Find the module-function (e.g. *IlluminationCorrectionPelkmans.m*).
2. Put two breakpoints:
	* at very start of the function:
		* before *handles* can be modified
		* the first *drawnow* statement is usually good
	* at very end of the function:
		* after *handles* is last modified
		* avoid sub-functions. Choose the final *end* statement of the main function.
3. Run *CellProfiler* and wait until the first break-point is reached.
	* At first break-point, execute
		> save( 'handles_in.mat', 'handles' );
 
	* Continue (press F5 as shortcut), and at second break-point, execute

		> save( 'handles_out.mat', 'handles' );

4. Look in your current working directory (type *pwd()* if unsure) for the files *handles\_in.mat* and *handles\_out.mat*. Check for expected file timestamps.

5. Copy these files into a new directory in:
	
	> **camelot-share-3**/Data/Datasets/GitHubTests/CellProfilerPelkmans/State

	Folders are grouped:
	
	1. The purpose of the test (e.g. a module name)
	2. Then by *A, B, C, etc.* referring to the associated input data
	3. Then by a number in the filename *1\_handles\_in.mat*, as there might be several cases. Please preserve the filename pattern.

6. Each of these folder has a file *pathOld.txt* that records the images-directory used when the snapshops were made. This will be recorded within the *handles* data structures. It is used for search-and-replace operations during the unit tests, to update the paths to a new location.