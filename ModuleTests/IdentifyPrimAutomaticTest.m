classdef IdentifyPrimAutomaticTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
            dirTestsSub = 'State/IdentifyPrimAutomatic';
            handlesIn = getTestHandles(dirTestsSub,'C/handles_in');
            handlesOut = getTestHandles(dirTestsSub,'C/handles_out');
            
            testResult = IdentifyPrimAutomatic( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end
    end
end