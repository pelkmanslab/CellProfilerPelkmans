function [ handles ] = getTestHandlesCopiedFolder( subFolderPath, handlesName, pathOld, pathNew, testCase )
% creates a temporary directory. copies test files to it, before rewriting
% paths
    tempFolderFixture = testCase.getSharedTestFixtures('matlab.unittest.fixtures.TemporaryFolderFixture');
    temporaryFolder = tempFolderFixture.Folder;

    % Copies all files from our test-folder (pathNew) to the temporary
    % folder
    copyfile(pathNew,temporaryFolder)
    
    % We rewrite paths to point to the copied folder
    handles = getRewrittenTestHandles( subFolderPath, handlesName, pathOld, temporaryFolder );
end
