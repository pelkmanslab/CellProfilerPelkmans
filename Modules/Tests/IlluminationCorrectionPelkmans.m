classdef IlluminationCorrectionPelkmans < matlab.unittest.TestCase
    
    properties
        A1_HandlesIn
        A1_HandlesOut
    end
 
    methods(TestMethodSetup)
        function createHandles(testCase)
            dirTestsSub = 'State/IlluminationCorrectionZScoring';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 3';
            pathNew = absolutePathFromRelative('Input/A');
            
            testCase.A1_HandlesIn = getRewrittenTestHandles(dirTestsSub,'A/1_handles_in', pathOld, pathNew);
            testCase.A1_HandlesOut = getRewrittenTestHandles(dirTestsSub,'A/1_handles_out', pathOld, pathNew);
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
            testResult = IlluminationCorrectionPelkmans( testCase.A1_HandlesIn );
            testCase.verifyTrue( isequal(testResult, testCase.A1_HandlesOut) );
        end
    end

end