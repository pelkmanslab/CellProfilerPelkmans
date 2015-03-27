function handles = SaveSegmentedCells_MPcycle(handles)

% Help for the SaveSegmentedCells module:
% Category: Other
%
% SHORT DESCRIPTION:
% Identifies objects (e.g. cell edges) using "seed" objects identified by
% an Identify Primary module (e.g. nuclei).
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Which objects do you want to save?
%choiceVAR01 = Do not use
%infotypeVAR01 = objectgroup
ObjectNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 =
%choiceVAR02 = Do not use
%infotypeVAR02 = objectgroup
ObjectNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 =
%choiceVAR03 = Do not use
%infotypeVAR03 = objectgroup
ObjectNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 =
%choiceVAR04 = Do not use
%infotypeVAR04 = objectgroup
ObjectNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Where do you want to save the segmentation images?
%choiceVAR05 = ReferenceCycle
%choiceVAR05 = CurrentCycle
OutputFolder = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu


%%%VariableRevisionNumber = 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    drawnow
    
    %%% retrieve shift descriptors from handles or ALIGNCYCLES directory
    if strcmp(getlastdir(handles.Current.DefaultOutputDirectory),'BATCH')
        strAlignCyclesDir = strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'],[filesep,'ALIGNCYCLES']);
    end
    if isfield(handles,'shiftDescriptor')
        shift = handles.shiftDescriptor;
        useShift = true;
    elseif isdir(strAlignCyclesDir)
        jsonDescriptorFile = fullfile(strAlignCyclesDir,'shiftDescriptor.json');
        shift = loadjson(jsonDescriptorFile);
        useShift = true;
    else
        useShift = false;
    end
    
    %%% get index of current image
    strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
    
    %%% get cycle number
    if ~isempty(strfind(handles.Current.DefaultImageDirectory,'TIFF'))
        iCycle = str2double(cell2mat(flatten(regexp(handles.Current.DefaultImageDirectory,'.*?([0-9]+).TIFF','tokens'))));
    else
        error(['Image processing was canceled in the ', ModuleName, ' module because the cycle number could not be determined correctly.'])
    end
    
    
    Objects = ObjectNameList(~strcmp(ObjectNameList,'Do not use'));
    
    for iObject = 1:length(Objects)
        
        strObjectName = char(Objects(iObject));
        
        fieldname = ['Segmented',strObjectName];
        
        matDotIndices = strfind(strOrigImageName,'.');
        % new CP apparently removes file extensions from image names
        if ~isempty(matDotIndices)
            strOrigImageName = strOrigImageName(1,1:matDotIndices(end)-1);
        end
        
        % if analyzing subdirectories, do not take along the subdirectory
        % structure... but replaces file separators with underscores to
        % preserver the additional information of the path.
        if ~isempty(strfind(strOrigImageName,filesep))
            strOrigImageName = strrep(strOrigImageName,filesep,'_');
        end

        if useShift  
            if strcmp(OutputFolder, 'ReferenceCycle')
                % built absolute SEGMENTATION path from relative path stored in handles
                strSegmentationDir = [strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], filesep), shift.SegmentationDirectory];
            else
                strSegmentationDir = strrep(handles.Current.DefaultOutputDirectory, 'BATCH', 'SEGMENTATION_ALIGNED');
            end
            
            [~, strMicroscopeType] = check_image_position(strOrigImageName);
            if strcmp(strMicroscopeType, 'CV7K')
                strExpression = '.+(_\w{1}\d{2}_)';
            elseif strcmp(strMicroscopeType, 'Visi')
                strExpression = '^.+(_s\d{4}_r\d{2}_c\d{2}_)';
            else
                error('Microscope type not implemented')
            end
            % built full segmentation image filename from filename trunk stored in handles
            if ~any(strcmp(strObjectName,{'Cells','Nuclei','Cytoplasm'}))
                % if object other than 'Cells', 'Nuclei' or 'Cytoplasm' append
                % filename with cycle number
                strSegmentationFileName = [regexprep(strOrigImageName, strExpression, sprintf('%s$1',shift.SegmentationFileNameTrunk)), ...
                                                    '_Segmented', strObjectName, sprintf('Cycle%d',iCycle), '.png'];
            else
                strSegmentationFileName = [regexprep(strOrigImageName, strExpression, sprintf('%s$1',shift.SegmentationFileNameTrunk)), ...
                                                    '_Segmented', strObjectName, '.png'];             
            end

            % determine output directory, if BATCH, create SEGMENTATION directory,
            % otherwise, use the default output directory to dump the segmentation
            % images.
            strOutputDir = strSegmentationDir;
            if strcmp(getlastdir(strOutputDir),'SEGMENTATION') || strcmp(getlastdir(strOutputDir),'SEGMENTATION_ALIGNED')
                if ~fileattrib(strOutputDir)
                    fprintf('%s: creating default output directory SEGMENTATION in %s',mfilename, getbasedir(strOutputDir))
                    mkdir(strOutputDir)
                end
            end

        else
            % built absolute SEGMENTATION path
            strSegmentationDir = strrep(handles.Current.DefaultOutputDirectory, 'BATCH', 'SEGMENTATION');

            % built full segmentation image filename
            if ~any(strcmp(strObjectName,{'Cells','Nuclei','Cytoplasm'}))
                % if object other than 'Cells', 'Nuclei' or 'Cytoplasm' append
                % filename with cycle number
                strSegmentationFileName = [strOrigImageName,'_Segmented',strObjectName,sprintf('Cycle%d',iCycle),'.png'];
            else
                strSegmentationFileName = [strOrigImageName,'_Segmented',strObjectName,'.png'];
            end

            strOutputDir = strSegmentationDir;
            
            if ~fileattrib(strOutputDir)
                fprintf('%s: creating default output directory SEGMENTATION in %s',mfilename,getbasedir(strOutputDir))
                mkdir(strOutputDir)
            end

        end
        
        strFullFileName = fullfile(strOutputDir, strSegmentationFileName);
        LabelMatrixImage = CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale');
        imwrite(uint16(LabelMatrixImage),strFullFileName,'png');
        
    end
    
end
