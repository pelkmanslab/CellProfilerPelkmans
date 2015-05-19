classdef MeasureObjectNeighborsPelkmansTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
            dirTestsSub = 'State/MeasureObjectNeighborsPelkmans';
            handlesIn = getTestHandles(dirTestsSub,'C/handles_in');
            handlesOut = getTestHandles(dirTestsSub,'C/handles_out');
            
            testResult = MeasureObjectNeighborsPelkmans( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end
    end
end