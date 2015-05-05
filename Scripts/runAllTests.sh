#!/bin/sh
#
# Bootstraps the doTesting.m matlab script
#
# Expects to be run from the top-level folder of CellProfilerPelkmans
#
# The Matlab script, however, will run in the folder Scripts/
#
filepath='./Scripts/runAllTests.m'
matlab -nodesktop -nosplash -nodisplay -r "run $filepath ; quit;" 