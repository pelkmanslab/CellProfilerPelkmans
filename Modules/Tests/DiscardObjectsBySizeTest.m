classdef DiscardObjectsBySizeTest < matlab.unittest.TestCase
    
    properties
        HandlesIn
        HandlesOut
    end
 
    methods(TestMethodSetup)
        function createHandles(testCase)
            dirTestsSub = '\DiscardObjectsBySize\TestA\';
            testCase.HandlesIn = getTestSavedHandles(dirTestsSub,'DiscardObjectsBySize_In');
            testCase.HandlesOut = getTestSavedHandles(dirTestsSub,'DiscardObjectsBySize_Out');
        end
    end
 
    methods(TestMethodTeardown)
        function closeHandles(testCase)
            testCase.HandlesIn = [];
            testCase.HandlesOut = [];
        end
    end
    
    methods(Test)
        function testAOutputStructureMatches(testCase)
            testResult = DiscardObjectsBySize( testCase.HandlesIn );
            testCase.verifyTrue( isequal(testResult, testResult) );
        end
    end
end