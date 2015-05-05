import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.TAPPlugin
import matlab.unittest.plugins.ToFile
%
% Script to execute matlab tests
%
% See:
% http://ch.mathworks.com/help/matlab/ref/matlab.unittest.plugins.tapplugin-class.html
%
% Expects to be run from the Tests/ directory
%

topLevelMatlabFolder = '../Modules/Tests';
outFilePathRelative = '../testsOutput.tap';

% All tests from the main library folder (recursively)
suites = matlab.unittest.TestSuite.fromFolder(topLevelMatlabFolder, 'IncludingSubfolders', true);

runner = TestRunner.withTextOutput;



% NOTE: It is critically important that an absolute file path is passed to
%   TAPPlugin.producingOriginalFormat as otherwise test entries become
%   missing in the output TAP file
outFileResolved = fullfile(pwd(),outFilePathRelative);

% We delete the existing tapFile as TAPPlugin.producingOriginalFormat
%   appends to any existing files, and we need it to be only a s
%    single test case for Jenkins.
delete(tapFile);
plugin = TAPPlugin.producingOriginalFormat(ToFile(outFileResolved));

runner.addPlugin(plugin)
result = runner.run(suites);

% Display TAP file
disp(fileread(outFileResolved))