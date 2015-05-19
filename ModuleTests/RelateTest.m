classdef RelateTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchC(testCase)
            
%             dirTestsSub = 'State/Relate';
%             handlesIn = getTestHandles(dirTestsSub,'C/handles_in');
%             handlesOut = getTestHandles(dirTestsSub,'C/handles_out');
%             
%             testResult = Relate( handlesIn );
%             testCase.verifyTrue( isequal(testResult, handlesOut) );
            
            

            
            
            dirTestsSub = 'State/Relate';
            
            pathOld = 'S:\Data\Users\Owen\TestDatasets\Dataset 5';
            pathNew = absolutePathFromRelative('Input/C');
            
            handlesIn = getRewrittenTestHandles(dirTestsSub,'C/handles_in', pathOld, pathNew);
            handlesOut = getRewrittenTestHandles(dirTestsSub,'C/handles_out', pathOld, pathNew);
            
            testResult = IdentifySpots2D( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );            
        end
    end
end