function [ handles ] = getTestHandlesRewrite( subFolderPath, handlesName, pathOld, pathNew )
% getTestSavedObject loads handles that have been serialized, and rewrites
% all occurences of pathOld to pathNew

    handles = getTestHandlesUnchanged( subFolderPath, handlesName );
    handles = struct_string_replace(  handles, pathOld, pathNew );
    handles = struct_string_replace(  handles, '\', '/' );

end

