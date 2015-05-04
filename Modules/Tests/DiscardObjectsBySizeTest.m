classdef DiscardObjectsBySizeTest < matlab.unittest.TestCase
    
    properties
        HandlesIn
        HandlesOut
    end
 
    methods(TestMethodSetup)
        function createHandles(testCase)
            
            dirTestsRoot = 'S:\Data\Users\Owen\CellProfilerPelkmansTestData';
            dirTestsSub = '\DiscardObjectsBySize\TestA\';
            
            load( fullfile(dirTestsRoot,dirTestsSub,'DiscardObjectsBySize_In.mat'),'handles');
            testCase.HandlesIn = handles;
            
            load( fullfile(dirTestsRoot,dirTestsSub,'DiscardObjectsBySize_Out.mat'),'handles');
            testCase.HandlesOut = handles;
        end
    end
 
    methods(TestMethodTeardown)
        function closeHandles(testCase)
            %close(testCase.TestFigure)
            testCase.HandlesIn = [];
            testCase.HandlesOut = [];
        end
    end
    
    methods(Test)
        function testAOutputStructureMatches(testCase)
            testResult = DiscardObjectsBySize( testCase.HandlesIn );
            testCase.verifyTrue( isequal(testResult, testResult) );
        end
    end
end