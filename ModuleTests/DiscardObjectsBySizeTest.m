classdef DiscardObjectsBySizeTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
            dirTestsSub = 'State/DiscardObjectsBySize';
            handlesIn = testHandlesUnchanged(dirTestsSub,'C/handles_in');
            handlesOut = testHandlesUnchanged(dirTestsSub,'C/handles_out');
            
            testResult = DiscardObjectsBySize( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end
    end
end