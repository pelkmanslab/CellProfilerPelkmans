import matlab.unittest.TestSuite
%
% Script to execute ALL matlab tests in the project
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
topLevelMatlabFolder = '../ModuleTests';
outFilePathRelative = '../testsOutput.tap';

% Adds the CellProfilerPelkmans base dir, and all sub-folders to the current path
mfilepath=fileparts(which(mfilename));
projectBaseDir=fullfile(mfilepath,'../');
addpath(genpath(projectBaseDir));

checkTestInit();

runner = createTestRunner(outFilePathRelative);

% All tests from the main library folder (recursively)
suites = matlab.unittest.TestSuite.fromFolder(topLevelMatlabFolder, 'IncludingSubfolders', true);
result = runner.run(suites);

% Display TAP file
disp(fileread(outFilePathRelative))