classdef DiscardObjectsBySizeTest < matlab.unittest.TestCase
    
    properties
        A1_HandlesIn
        A1_HandlesOut
    end
 
    methods(TestMethodSetup)
        function createHandles(testCase)
            dirTestsSub = 'State/DiscardObjectsBySize';
            testCase.A1_HandlesIn = getTestHandles(dirTestsSub,'C/1_handles_in');
            testCase.A1_HandlesOut = getTestHandles(dirTestsSub,'C/1_handles_out');
        end
    end
 
    methods(TestMethodTeardown)
        function closeHandles(testCase)
            testCase.A1_HandlesIn = [];
            testCase.A1_HandlesOut = [];
        end
    end
    
    methods(Test)
        function handlesMatchA(testCase)
            testResult = DiscardObjectsBySize( testCase.A1_HandlesIn );
            testCase.verifyTrue( isequal(testResult, testCase.A1_HandlesOut) );
        end
    end
end