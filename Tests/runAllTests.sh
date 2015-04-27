#!/bin/sh
#
# Bootstraps the doTesting.m matlab script
#
# Expects to be run from the top-level folder of CellProfilerPelkmans
#
$filename = 'runAllTests.m'
matlab -nodesktop -nosplash -nodisplay -r "run ./$filename ; quit;" 