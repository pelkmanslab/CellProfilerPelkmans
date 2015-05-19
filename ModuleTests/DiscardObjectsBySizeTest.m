classdef DiscardObjectsBySizeTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
            dirTestsSub = 'State/DiscardObjectsBySize';
            handlesIn = getTestHandles(dirTestsSub,'C/1_handles_in');
            handlesOut = getTestHandles(dirTestsSub,'C/1_handles_out');
            
            testResult = DiscardObjectsBySize( testCase.handlesIn );
            testCase.verifyTrue( isequal(testResult, testCase.handlesOut) );
        end
    end
end