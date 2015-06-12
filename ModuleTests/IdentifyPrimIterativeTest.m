classdef IdentifyPrimIterativeTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchD(testCase)
            
            dirTestsSub = 'State/IdentifyPrimIterative';
            handlesIn = testHandlesUnchanged(dirTestsSub,'D/handles_in');
            handlesOut = testHandlesUnchanged(dirTestsSub,'D/handles_out');
            
            testResult = IdentifyPrimIterative( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end

    end
end
