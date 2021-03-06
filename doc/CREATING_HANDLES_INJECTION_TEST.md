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
	2. Then by *A, B, C, etc.* referring to the associated input data (From Step 1)
	3. Then by *handles\_in.mat* or *handles\_out.mat*.

6. Store the *recording-path* in a file *pathOld.txt* that records the images-directory used when the snapshots were made. This is:
	* **NB: This step is important.**
	* The exact path used when the state was recorded e.g.
		> S:\Data\Users\Owen\TestDatasets\Dataset 5
		
	* It is the folder where the pipeline resides (e.g. parent folder of *TIFF/*)
	* It will vary by operating-system.
	* Do not include a trailing /
	
## Step 3: Write a verification test

1. Add a new test in *ModuleTests* following the *ModuleName*Test.m filename pattern.
2. Copy an existing-test and update the functions and paths accordingly. Choose a suitable template test, as explained in the subsequent section.
3. Check that the test works locally, by running *runSpecificTest.m* or *runAllTests.m* with *Scripts/* as your current working directory.
4. Commit and sync with GitHub. Jenkins will run the test remotely. You will get an email if any tests fail.

## Choosing a template test

Modules modify the *handles* structure-array. But they vary in their file I/O, with three cases:

1. **Some don't do any file I/O.** Use *getTestHandlesUnchanged* function. See  *DiscardObjectsBySizeTest.m*
2. **Some read files.** Use *getTestHandlesRewrite* function. See  *IlluminationCorrectionPelkmansTest.m*
3. **Some read/write files.** Use *getTestHandlesCopy* function. See  *SaveSegmentedObjectsTest.m*

Choose the correct case based upon knowledge of the module. Or trial and error, in the above order.

Path-rewriting needs to occur in Case 2 and 3.

We consider the test-data as read-only; Jenkins explicitly mounts test-data as read-only. Accordingly for Case 3, a local copy of the test-data is first-made, and deleted afterwards. These tests also need a *TemporaryFolderFixture* in Matlab.