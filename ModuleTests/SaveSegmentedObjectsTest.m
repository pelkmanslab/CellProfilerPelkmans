classdef SaveSegmentedObjectsTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
      
            dirTestsSub = 'State/SaveSegmentedObjects';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 5';
            pathNew = absolutePathFromRelative('Input/C');
            
            handlesIn = getRewrittenTestHandles(dirTestsSub,'C/handles_in', pathOld, pathNew);
            handlesOut = getRewrittenTestHandles(dirTestsSub,'C/handles_out', pathOld, pathNew);
            
            testResult = SaveSegmentedObjects( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );            
        end
    end
end