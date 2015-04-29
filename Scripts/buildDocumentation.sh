#!/bin/sh
#
# Bootstraps the buildDocumentation.m matlab script
#
# Expects to be run from the top-level folder of iBRAINShared
#
#
filepath='./Scripts/buildDocumentation.m'
matlab -nodesktop -nosplash -nodisplay -r "run $filepath ; quit;" 


