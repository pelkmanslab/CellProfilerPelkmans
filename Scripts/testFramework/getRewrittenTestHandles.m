function [ handles ] = getRewrittenTestHandles( subFolderPath, handlesName, pathOld, pathNew )
% getTestSavedObject loads handles that have been serialized, and rewrites
% all occurences of pathOld to pathNew

    handles = getTestHandles( subFolderPath, handlesName );
    handles = struct_string_replace(  handles, pathOld, pathNew );

end

