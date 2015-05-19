classdef IlluminationCorrectionPelkmansTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchA(testCase)
            
            dirTestsSub = 'State/IlluminationCorrectionZScoring';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 3';
            pathNew = absolutePathFromRelative('Input/A');
            
            handlesIn = testHandlesRewrite(dirTestsSub,'A/handles_in', pathOld, pathNew);
            handlesOut = testHandlesRewrite(dirTestsSub,'A/handles_out', pathOld, pathNew);
            
            testResult = IlluminationCorrectionPelkmans( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end
    end

end