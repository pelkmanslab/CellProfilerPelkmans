classdef (SharedTestFixtures={matlab.unittest.fixtures.TemporaryFolderFixture}) ...
    SaveSegmentedObjectsTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
      
            dirTestsSub = 'State/SaveSegmentedObjects';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 5';
            pathNew = absolutePathFromRelative('Input/C');
            
            handlesIn = testHandlesCopy(dirTestsSub,'C/handles_in', pathOld, pathNew, testCase);
            handlesOut = testHandlesCopy(dirTestsSub,'C/handles_out', pathOld, pathNew, testCase);
            
            testResult = SaveSegmentedObjects( handlesIn );
            testHandlesVerifyEqual(testCase, testResult, handlesOut);
        end
    end
end