# Creating a Handles-Injection Test

A handles-injection test is a particular type of check performed on CellProfiler modules. 

CellProfiler modules accept a structure-array *handles* as an argument. Some processing is performed which alters *handles* and it is returned.
 
This type of test verifies that the behaviour of a module remains consistent by recording and comparing examples of *handles*.

There are three steps necessary:

1. **Extract test data**. Find an example pipeline which is known to function correctly (for the module of interest). Select a small subset of images.
2. **Record handles-state**. Run this pipeline locally on your PC. Record *handles* when entering and leaving the module.
3. **Write a verification test**. Verify that the module produces the same output in future.

Basically, the check performed is:

> module\_function( *recorded\_handles\_in\_state* ) = *recorded\_handles\_out\_state*

In practice, the some rewriting of paths needs to occur, as the filesystem paths stored in *handles* will differ between environments:

* **recording environment** - where in-state and out-state are recorded (one particular personal computer)
* **testing environment** - for testing in the future (on Jenkins server, varying personal computers) 

So the **recording-path** is saved with the in-state and out-state. Paths are then rewritten as needed with a *search-and-replace* operation.

## Step 1: Extract test data 


## Step 2: Record handles-state

1. Find the module-function (e.g. IlluminationCorrectionPelkmans.m).
2. Put breakpoints
	* at very start of the function (before *handles* can be modified).
	* at very end of the function (after *handles* is last modified).
3. Run *CellProfiler* and wait until the first break-point is reached.
	* At first break-point, execute
		> save( 'handles_in.mat', 'handles' );
 
	* Resume, and at second break-point, execute

		> save( 'handles_out.mat', 'handles' );

4.