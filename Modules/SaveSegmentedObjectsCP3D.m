function handles = SaveSegmentedObjectsCP3D(handles)

% Help for the SaveSegmentedObjectsCP3D module:
% Category: Other
%
% SHORT DESCRIPTION:
% Saves the Segmentation of CP3D Objects
%
%   Compatibility: largely follows code of SaveSegmentedCells, except that
%   it checks if there is .Format within the object (which is false for
%   cp1, but true when following stanadard of cp3d) in that case the output
%   will be the structure of the object (which dependent upon the format
%   choosen may contain the linear indices of individual objects)
%   FILENAMES are derived from the current Measuremnt.Image which, for 3D
%   has a sole filename per stack.
%
% *************************************************************************
%
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});


%textVAR02 = 


%%%VariableRevisionNumber = 145

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

fieldname = ['Segmented',strObjectName];

% check, if CP3D Segmentation Information
if isfield(handles.Pipeline.(fieldname),'Format') == true
    bnHasCP3DSegmentation = true;
else
    bnHasCP3DSegmentation = false;
end

% create new image name
% compatibility with both old and new CellProfiler measurement storing
% schemas.
if isfield(handles.Measurements.Image, 'FileNames')
    strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
else
    cellstrImageFieldNames = fieldnames(handles.Measurements.Image);
    strImageFieldName = cellstrImageFieldNames{find(strncmp(cellstrImageFieldNames,'FileName_',9),1,'first')};
    strOrigImageName = char(handles.Measurements.Image.(strImageFieldName){handles.Current.SetBeingAnalyzed});
end

matDotIndices=strfind(strOrigImageName,'.');
% new CP apparently removes file extensions from image names
if ~isempty(matDotIndices)
    strOrigImageName = strOrigImageName(1,1:matDotIndices(end)-1);
end

switch bnHasCP3DSegmentation
    case false
        strTmpFileName = [strOrigImageName,'_',fieldname,'.png'];
    case true
        outName = ['SegmentedCP3D',strObjectName];
        strTmpFileName = [strOrigImageName,'_',outName,'.mat'];
end


% if analyzing subdirectories, do not take along the subdirectory
% structure... but replaces file separators with underscores to
% preserver the additional information of the path.
if ~isempty(strfind(strOrigImageName,filesep))
    strTmpFileName = strrep(strTmpFileName,filesep,'_');
end

% determine output directory, if BATCH, create SEGMENTATION directory,
% otherwise, use the default output directory to dump the segmentation
% images.
strOutputDir = handles.Current.DefaultOutputDirectory;
if strcmp(getlastdir(strOutputDir),'BATCH')
    strOutputDir = strrep(strOutputDir, [filesep,'BATCH'],[filesep,'SEGMENTATION']);
    if ~fileattrib(strOutputDir)
        disp(sprintf('%s: creating default output directory SEGMENTATION in %s',mfilename,getbasedir(strOutputDir)))
        mkdir(strOutputDir)
    end
end

strFullFileName = fullfile(strOutputDir, strTmpFileName);

switch bnHasCP3DSegmentation
    case false
        LabelMatrixImage = CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale');
        imwrite(uint16(LabelMatrixImage),strFullFileName,'png');
    case true
        LabelStruct =  handles.Pipeline.(fieldname);    %LabelStruct will be saved as a variable
        save(strFullFileName, 'LabelStruct','-mat');
end

end
