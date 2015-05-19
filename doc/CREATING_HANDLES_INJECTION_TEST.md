# Creating a Handles-Injection Test

A handles-injection test is a particular type of check performed on CellProfiler modules. 

CellProfiler modules accept a structure-array *handles* as an argument. Some processing is performed which alters *handles* and it is returned.
 
This type of test verifies that the behaviour of a module remains consistent.  There are three steps necessary:

1. **Extract test data**. Find an example pipeline which is known to function correctly (for the module of interest). Select a small subset of images.
2. **Record handles-state**. Run this pipeline locally on your PC. Record *handles* when entering and leaving the module.
3. **Write a verification test**. Verify that the module produces the same output in future.

Basically, the check performed is:

> module\_function( *recorded\_handles\_in* ) = *recorded\_handles\_out*

## Step 1: Extract test data 