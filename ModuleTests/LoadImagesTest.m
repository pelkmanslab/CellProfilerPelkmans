classdef LoadImagesTest < matlab.unittest.TestCase
    
    properties
        A1_HandlesIn
        A1_HandlesOut
    end
 
    methods(TestMethodSetup)
        function createHandles(testCase)
            dirTestsSub = 'State/LoadImages';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 3\TIFF';
            pathNew = absolutePathFromRelative('Input/A/TIFF');
            
            testCase.A1_HandlesIn = testHandlesRewrite(dirTestsSub,'A/handles_in', pathOld, pathNew);
            testCase.A1_HandlesOut = testHandlesRewrite(dirTestsSub,'A/handles_out', pathOld, pathNew);

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
            testHandlesVerifyEqual(testCase, testResult, testCase.A1_HandlesOut);
        end
    end

end