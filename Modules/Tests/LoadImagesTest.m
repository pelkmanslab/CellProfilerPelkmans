classdef LoadImagesTest < matlab.unittest.TestCase
    
    properties
        A1_HandlesIn
        A1_HandlesOut
    end
 
    methods(TestMethodSetup)
        function createHandles(testCase)
            dirTestsSub = 'LoadImages';
            testCase.A1_HandlesIn = getTestSavedHandles(dirTestsSub,'A1_handles_in');
            testCase.A1_HandlesOut = getTestSavedHandles(dirTestsSub,'A1_handles_out');
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
            testResult = LoadImages( testCase.A1_HandlesIn );
            testCase.verifyTrue( isequal(testResult, testCase.A1_HandlesOut) );
        end
    end
end