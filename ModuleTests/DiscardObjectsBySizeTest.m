classdef DiscardObjectsBySizeTest < matlab.unittest.TestCase
    
    properties
        C1_HandlesIn
        C1_HandlesOut
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
        function handlesMatchC(testCase)
            testResult = DiscardObjectsBySize( testCase.C1_HandlesIn );
            testCase.verifyTrue( isequal(testResult, testCase.C1_HandlesOut) );
        end
    end
end