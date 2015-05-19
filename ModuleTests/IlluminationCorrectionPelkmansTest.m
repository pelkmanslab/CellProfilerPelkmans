classdef IlluminationCorrectionPelkmansTest < matlab.unittest.TestCase
    
    properties
        handlesIn
        handlesOut
    end
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchA(testCase)
            
            dirTestsSub = 'State/IlluminationCorrectionZScoring';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 3';
            pathNew = absolutePathFromRelative('Input/A');
            
            handlesIn = getRewrittenTestHandles(dirTestsSub,'A/handles_in', pathOld, pathNew);
            handlesOut = getRewrittenTestHandles(dirTestsSub,'A/handles_out', pathOld, pathNew);
            
            testResult = IlluminationCorrectionPelkmans( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end
    end

end