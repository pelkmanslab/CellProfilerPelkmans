classdef IdentifySpots2DTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
%             dirTestsSub = 'State/IdentifySpots2D';
%             handlesIn = getTestHandles(dirTestsSub,'C/handles_in');
%             handlesOut = getTestHandles(dirTestsSub,'C/handles_out');
%             
%             testResult = IdentifySpots2D( handlesIn );
%             testCase.verifyTrue( isequal(testResult, handlesOut) );
            
            
            dirTestsSub = 'State/IdentifySpots2D';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 5';
            pathNew = absolutePathFromRelative('Input/C');
            
            handlesIn = testHandlesRewrite(dirTestsSub,'C/handles_in', pathOld, pathNew);
            handlesOut = testHandlesRewrite(dirTestsSub,'C/handles_out', pathOld, pathNew);
            
            testResult = IdentifySpots2D( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end
    end
end