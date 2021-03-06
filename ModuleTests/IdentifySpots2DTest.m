classdef IdentifySpots2DTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
           
            dirTestsSub = 'State/IdentifySpots2D';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 5';
            pathNew = absolutePathFromRelative('Input/C');
            
            handlesIn = testHandlesRewrite(dirTestsSub,'C/handles_in', pathOld, pathNew);
            handlesOut = testHandlesRewrite(dirTestsSub,'C/handles_out', pathOld, pathNew);
            
            testResult = IdentifySpots2D( handlesIn );
            testHandlesVerifyEqual( testCase, testResult, handlesOut );
        end
    end
end