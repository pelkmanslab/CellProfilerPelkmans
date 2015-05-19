function [ handlesOut ] = getTestHandlesUnchanged( subFolderPath, handlesName )
% getTestSavedObject loads handles that have been serialized
%
%   loads a file from the test directory, using
%   baseFolder/subFolderPath/fileName.mat
    
   fileName = strcat(handlesName,'.mat');
   fullFilePath = fullfile( absolutePathFromRelative(subFolderPath),fileName);
   load( fullFilePath,'handles');
   handlesOut = handles;
end

