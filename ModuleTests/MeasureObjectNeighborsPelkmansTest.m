classdef MeasureObjectNeighborsPelkmansTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
            dirTestsSub = 'State/MeasureObjectNeighborsPelkmans';
            handlesIn = testHandlesUnchanged(dirTestsSub,'C/handles_in');
            handlesOut = testHandlesUnchanged(dirTestsSub,'C/handles_out');
            
            testResult = MeasureObjectNeighborsPelkmans( handlesIn );
            testHandlesVerifyEqual(testCase, testResult, handlesOut);
        end
    end
end