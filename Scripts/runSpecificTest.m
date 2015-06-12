import matlab.unittest.TestSuite
%
% Script to execute a specific matlab test from a XXXXXTest.m
%
% See:
% http://ch.mathworks.com/help/matlab/ref/matlab.unittest.plugins.tapplugin-class.html
%
% Expects to be run from the Scripts/ directory
%
% An environmental variable CELLPROFILER_PELKMANS_TESTS_DIR must be set
%   with location of test data
%
% Results are outputted in TAP format to outFilePathRelative
%
testFilePathRelative = '../ModuleTests/IdentifyPrimIterativeTest.m';
outFilePathRelative = '../testsOutput.tap';

% Adds the CellProfilerPelkmans base dir, and all sub-folders to the current path
mfilepath=fileparts(which(mfilename));
projectBaseDir=fullfile(mfilepath,'../');
addpath(genpath(projectBaseDir));

checkTestInit();

runner = createTestRunner(outFilePathRelative);

suite = matlab.unittest.TestSuite.fromFile(testFilePathRelative);
result = runner.run(suite);

% Display TAP file
disp(fileread(outFilePathRelative))
