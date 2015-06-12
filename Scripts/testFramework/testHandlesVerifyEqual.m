function testHandlesVerifyEqual( testCase, testResult, handlesOut )
%TESTHANDLESVERIFYEQUAL Verifies that two handles objects are equal
%   Detailed explanation goes here

    match = isequal(testResult, handlesOut);
    
    if match==0 
        disp('\n\nNot matching. Starting a diff of the two objects:');        
        structcmp(testResult, handlesOut, 'Report', 'on')
    end
    
    testCase.verifyTrue( match );
end

