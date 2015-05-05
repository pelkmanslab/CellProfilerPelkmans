function [ handlesOut ] = getTestSavedHandles( subFolderPath, handlesName )
% getTestSavedObject loads handles that have been serialized
%
%   loads a file from the test directory, using
%   baseFolder/subFolderPath/fileName
   basePath = getenv('CELLPROFILER_PELKMANS_TESTS_DIR');
    
   fileName = strcat(handlesName,'.mat');
   fullFilePath = fullfile(basePath,subFolderPath,fileName);
   load( fullFilePath,'handles');
   handlesOut = handles;
end

