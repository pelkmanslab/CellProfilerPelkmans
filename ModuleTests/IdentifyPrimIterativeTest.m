classdef IdentifyPrimIterativeTest < matlab.unittest.TestCase
 
    methods(TestMethodSetup)
    end
 
    methods(TestMethodTeardown)
    end
    
    methods(Test)
        function handlesMatchD(testCase)
            
            dirTestsSub = 'State/IdentifyPrimIterative';

            pathOld = '/Users/mdh/jtui/expdata/150316-30min-PBS/150316-30min-PBS_03';
            pathNew = absolutePathFromRelative('Input/D');
            
            handlesIn = testHandlesRewrite(dirTestsSub,'D/handles_in', pathOld, pathNew);
            handlesOut = testHandlesRewrite(dirTestsSub,'D/handles_out', pathOld, pathNew);
        
            testResult = IdentifyPrimIterative( handlesIn );
            testCase.verifyTrue( isequal(testResult, handlesOut) );
        end

    end
end
