classdef IdentifyPrimAutomaticTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
            dirTestsSub = 'State/IdentifyPrimAutomatic';
            handlesIn = testHandlesUnchanged(dirTestsSub,'C/handles_in');
            handlesOut = testHandlesUnchanged(dirTestsSub,'C/handles_out');
            
            testResult = IdentifyPrimAutomatic( handlesIn );
            testHandlesVerifyEqual(testCase, testResult, handlesOut);
        end
    end
end