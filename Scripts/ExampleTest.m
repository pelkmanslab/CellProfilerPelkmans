classdef ExampleTest < matlab.unittest.TestCase
    methods(Test)
        function dummyTestOne(testCase)  % Test fails
            testCase.verifyEqual(5, 4, 'Testing 5==4')
        end
        function dummyTestTwo(testCase)  % Test passes
            testCase.verifyEqual(5, 5, 'Testing 5==5')
        end
        function dummyTestThree(testCase)
            % test code
        end
    end
end