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

suite   = TestSuite.fromClass(?ExampleTest);

runner = TestRunner.withTextOutput;

tapFile = '../testsOutput.tap';
plugin = TAPPlugin.producingOriginalFormat(ToFile(tapFile));

runner.addPlugin(plugin)
result = runner.run(suite);

disp(fileread(tapFile))