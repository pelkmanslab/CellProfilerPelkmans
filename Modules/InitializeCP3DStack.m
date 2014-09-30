function handles = InitializeCP3DStack(handles)
% Help for the InitializeCP3DStack module:
% Category: File Processing
%
%
% IMPORTANT NOTE:
% This module gets information about the metainformation from the file name
% of the images. These are microscope specific. Unless you use exactly the
% same microsocpe, you might have to adjust METRFROMIMAGENAME and there
% write your own custom code, which parses the file name and extracts the
% metainformation. You might have a look at matlab's help on regular
% expressions.
%
%
%
%
% SHORT DESCRIPTION:
% Allows you to initialize stacks based upon the metainformation contained
% in the filenames of individual images. This module does not directly load
% the images of the stacks into memory to allow a larger flexibility during
% the pipeline. Place LOADCP3DSTACK below this INITIALIZECP3DSTACK module
% to load the images into memory and UNLOADCP3DSTACK to remove them from
% memory. Manually controlling the memory allows to keep the memory
% requirment low, thus allowing parallelization on thousands of nodes with
% relatively low memory/node. 
%
% This module is basically a heavily modified version of the orignal
% LoadImages module of CP.
%   
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
%



%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = Type the text that one type of image has in common (for TEXT options):
%defaultVAR01 = C01.
TextToFind{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What do you want to call these images within CellProfiler?
%defaultVAR02 = StackBlue
%infotypeVAR02 = imagegroup indep
StackName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which Z positions do you want to load, 'all' will load?
%defaultVAR03 = all
ZPositions{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Type the text that one type of image has in common (for TEXT options):
%defaultVAR04 = /
TextToFind{2} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = What do you want to call these images within CellProfiler?
%defaultVAR05 = /
%infotypeVAR05 = imagegroup indep
StackName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Which Z positions do you want to load, 'all' will load?
%defaultVAR06 = /
ZPositions{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Type the text that one type of image has in common (for TEXT options):
%defaultVAR07 = /
TextToFind{3} = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = What do you want to call these images within CellProfiler?
%defaultVAR08 = /
%infotypeVAR08 = imagegroup indep
StackName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Which Z positions do you want to load, 'all' will load?
%defaultVAR09 = /
ZPositions{3} = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = Type the text that one type of image has in common (for TEXT
%options):
%defaultVAR10 = /
TextToFind{4} = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = What do you want to call these images within CellProfiler?
%defaultVAR11 = /
%infotypeVAR11 = imagegroup indep
StackName{4} = char(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = Which Z positions do you want to load, 'all' will load?
%defaultVAR12 = /
ZPositions{4} = char(handles.Settings.VariableValues{CurrentModuleNum,12});

%textVAR13 = Type the text that one type of image has in common (for TEXT options):
%defaultVAR13 = /
TextToFind{5} = char(handles.Settings.VariableValues{CurrentModuleNum,13});

%textVAR14 = What do you want to call these images within CellProfiler?
%defaultVAR14 = /
%infotypeVAR14 = imagegroup indep
StackName{5} = char(handles.Settings.VariableValues{CurrentModuleNum,14});

%textVAR15 = Which Z positions do you want to load, 'all' will load?
%defaultVAR15 = /
ZPositions{5} = char(handles.Settings.VariableValues{CurrentModuleNum,15});

%pathnametextVAR16 = Enter the path name to the folder where the images to be loaded are located. Type period (.) for default image folder.
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,16});


%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




drawnow

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%% Remove slashes entries with N/A or no filename from the input,
%%% i.e., store only valid entries
tmp1 = {};
tmp2 = {};
tmp3 = {};
for n = 1:length(StackName)
    if ~strcmp(TextToFind{n}, '/') && ~strcmp(StackName{n}, '/') && ~strcmp(ZPositions{n}, '/')
        tmp1{end+1} = TextToFind{n};
        tmp2{end+1} = StackName{n};
        tmp3{end+1} = ZPositions{n};
    end
end
TextToFind = tmp1;
StackName = tmp2;
ZPositions = tmp3;
clear tmp1; clear tmp2; clear tmp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  INITIALIZE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnUseRefNames = false;      % indicates whether image names initialized by an LoadImages module should be used as a reference. This is important to ensure that z-stacks and 2D images are grouped together correctly.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numGroups = length(TextToFind);

if SetBeingAnalyzed == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SYNCHRONIZE WITH LOADING MODULES WITHIN THE PIPELINE  %%%%%%
    % This will check, if any loading module is present in the pipeline. To
    % prevent incompatibilities, the initializeCP3D stack module has to be
    % downstream of the the LoadImage module. The first image of the load
    % image module will serve as a reference
    %
    %   ! ATTENTION ! This checking is very dependent upon the design of
    %   the CP LoadImages module. Several references are hardcoded.
    %   Similarly the check for a loading module works via naming the the
    %   names of the Moduls
    %
    
    f = cell2mat(cellfun(@(x) ~isempty(x), regexp(handles.Settings.ModuleNames,'LoadImages|LoadMoreImages|LoadEvenMoreImages)'), 'UniformOutput',0));
    if any(f)   % if any of the Loading Modules specified one line above, is present. Currently this excludes the laodindividualimage  module
        fIX = find(f,1,'first');  % obtain module number of first loading module
        if fIX>CurrentModuleNum     % let crash if initialize CP3D module occurs before load images module
            error(['Image processing was canceled in the ', ModuleName, ' module because it must not be placed upstream of a LoadImages module within the pipeline to ensure compatibility of 3D and 2D data.'])
        else
            names = fieldnames(handles.Pipeline);
            fi = cell2mat(cellfun(@(x) ~isempty(x), strfind(names,'Filename'),'UniformOutput',0));
            names = names(fi);
            
            strCandFilenameImages = ['Filename' handles.Settings.VariableValues{fIX,3}];      % field 3 corresponds to name of image !
            
            if any(ismember(names,strCandFilenameImages))
                strCandFileList = ['FileList' handles.Settings.VariableValues{fIX,3}];
                RefFileList = handles.Pipeline.(strCandFileList)';          % note the transposition: CP3D lists filenames of imagegroups(stacks) along rows, whereas CP along columns
                bnUseRefNames = true;
            else
                error(['Image processing was canceled in the ', ModuleName, ' module because it could not find reference filenames: Variable 3 of first LoadImages module must be set. '])
            end
        end
        % else bnUseRefNames will still be false
    end
    % Load metainformation of the reference images
    if bnUseRefNames == true
        [usubRCTSRef,subMaxRef,~] = obtainRCTS(RefFileList);
        if length(unique(usubRCTSRef))~=length(usubRCTSRef)
            error(['Image processing was canceled in the ', ModuleName, ' module because images previously loaded by Load Images are not unambiguous for a combination of row, column, site and timepoint. '])
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% START ANALYSING IMAGES FOR CP3D STACK  %%%%%%%%%%%%%%%%%%%%%
    
    %%%%% SET INPUT DIRECTORY %%%%%%
    if strncmp(Pathname,'.',1)
        if length(Pathname) == 1
            Pathname = handles.Current.DefaultImageDirectory;
        else
            Pathname = fullfile(handles.Current.DefaultImageDirectory,Pathname(2:end));
        end
    end
    SpecifiedPathname = Pathname;   % Check that input directory exists
    if ~exist(SpecifiedPathname,'dir')
        error(['Image processing was canceled in the ', ModuleName, ' module because the directory "',SpecifiedPathname,'" does not exist. Be sure that no spaces or unusual characters exist in your typed entry and that the pathname of the directory begins with /.'])
    end
    
    %%%%% PROCESS STACK-GROUPS %%%%%%
    for k=1: numGroups   % loop through all specified stack-groups
        
        %%% Load Filenames
        FileList = CPretrievemediafilenames(SpecifiedPathname,char(TextToFind{k}),'', 'Regular','Image')'; 
       
        %Follow logic of LoadImage module, But very SLOW!, does not even use fast CP modules, which replace standard matlab modules. Takes around 100 sec for a 60 wells. In  addition slow part of the function is repeated in this loop.
        
        %%% Checks whether any files are left.
        if isempty(FileList)
            error(['Image processing was canceled in the ', ModuleName, ' module because there are no image files with the text "', TextToFind{k}, '" in the chosen directory (or subdirectories, if you requested them to be analyzed as well).'])
        end
        
        %%% Extract Metainormation from these Filenames and Group according
        %%% to Combination of RowColumnTimepointSite , note that IDs of
        %%% intput (TextToFind) should mark different channels
        [subRCTS,subMax,intZPos] = obtainRCTS(FileList);
        
        %%% check if range of reference images (load image module) and of
        %%% stackimages is different
        if bnUseRefNames == true;
            if any(subMax ~= subMaxRef)
                error(['Image processing was canceled in the ', ModuleName, ' module because because images of the files with the text "', TextToFind{k}, '" are from different row or column or site or timepoint than the images loaded by a load images module.'])
            end
        end
        
        %%% define z-planes
        % In case that all Z positions should be loaded, set numerical
        % values of Z planes as span of smallest to largest otherwise check
        % if input consists of allowe characters and make vector with these values
        if strcmp(ZPositions{k},'all') == true
            ZPositions{k}=unique(intZPos); % if 'all' is specified, go through eall z planes in ascending order. make
        else
            % check if input only contains allowed characters
            [isSafe,ZPositions{k}]= inputVectorsForEvalCP3D(ZPositions{k},false);
            if isSafe == false
                error(['Image processing was canceled in the ', ModuleName, ' module because ' ZPositions{k} ' is not allowed input'])
            else
                % create vector
                ZPositions{k} = eval(ZPositions{k});
            end
        end
        % Make ZPositions as vector 1xn ; this will be important so that
        % column in metainformation will correspond to column of filenames
        if size(ZPositions{k},1)>size(ZPositions{k},2)
            ZPositions{k} = ZPositions{k}';
        end
        % Check if z planes are not specified multiple tiems
        if length(unique(ZPositions{k}))~=length(ZPositions{k})
            error(['Image processing was canceled in the ', ModuleName, ' module because some of ' ZPositions{k} ' of ' StackName{k} ' are present multiple times'])
        end
        
        %%% Filter data for requested z-planes %%%
        f = ismember(intZPos,ZPositions{k});
        subRCTS = subRCTS(f,:);
        intZPos = intZPos(f);
        FileList = FileList(f);
        
        % Check if all the requested z planes are present
        if any(~ismember(ZPositions{k},intZPos)) == true
            error(['Image processing was canceled in the ', ModuleName, ' module because not all ' ZPositions{k} ' Z planes could be found for ' StackName{k} '.'])
        end
        
        %%% Obtain IDs of images to group in one stack %%%
        [usubRCTS,uPos]= unique(subRCTS);       % obtain unique subRCTS
        [~, sortIX]=sort(uPos);                 % restore order of initial file-list > usubRCTSi
        usubRCTSi = usubRCTS(sortIX);
        
        % Check if correct amount of images
        ImgPerRCTS = histc(subRCTS,usubRCTS); % Note that edge vector must be monotonically non-decreasing (thus usubRCTS rather than usubRCTSi is used.
        if any(ImgPerRCTS ~= length(ZPositions{k})) == true
            error(['Image processing was canceled in the ', ModuleName, ' module because the number of files specified by "', TextToFind{k}, '" does not match the number of the requested Z-planes for each group of images.'])
        end
        clear usubRCTS;  % clear usubRCTS so that it will not be available later, making sure that only unsorted usubRCTSi will be used in later code.
        
        %%% SYNCHRONIZE WITH METADATA OF REFERENCE DATASETS %%%
        % this will ensure that images of one cycle will refer to the same
        % conbination of row/column/site/timepoint
        if bnUseRefNames == true
            if any(size(usubRCTSi) ~=  size(usubRCTSRef))   % size must be same
                error(['Image processing was canceled in the ', ModuleName, ' module because the number of row/column/site/time groups specified by "', TextToFind{k}, '" does not correspond to previously initialized (or loaded) groups.'])
            elseif any(usubRCTSi ~= usubRCTSRef) % if Reference and input usubRCTS are not identical, use reference
                if length(unique(usubRCTSi)) ~= length(unique(usubRCTSRef))
                    error(['Image processing was canceled in the ', ModuleName, ' module because the number of row/column/site/time groups specified by "', TextToFind{k}, '" does not correspond to previously initialized (or loaded) groups.'])
                elseif any(~ismember(usubRCTSi,usubRCTSRef))    % everhting must be mapped
                    error(['Image processing was canceled in the ', ModuleName, ' module because at least one of row/column/site/time groups specified by "', TextToFind{k}, '" does not correspond to any previously initialized (or loaded) group.'])
                end
            end
        else
            % if there is no preexisting reference, set first processed group of
            % this module as reference
            usubRCTSRef = usubRCTSi;
            subMaxRef = subMax;
            bnUseRefNames = true;
        end
        
        %%% INITIALIZE GROUPS OF Z-STACKS  %%%
        % Group Filenames, rows will indicate different z stacks and
        % columns different z-layers
        dimImageStacks = size(usubRCTSi,1);
        dimImagesPerStack = size(ZPositions{k},2);
        
        StackFileList = cell(dimImageStacks,dimImagesPerStack);
        for l=1:dimImageStacks
            for r = 1: dimImagesPerStack
                currentusubRCTS = usubRCTSi(l,1);
                f = (subRCTS == currentusubRCTS) & (intZPos==ZPositions{k}(1,r));
                StackFileList(l,r)=FileList(f,1);
            end
        end
        
        % Save information to handles
        fieldname = ['FileList', StackName{k}];
        handles.Pipeline.(fieldname) = StackFileList;
        fieldname = ['Pathname', StackName{k}];
        handles.Pipeline.(fieldname) = SpecifiedPathname;
        fieldname = ['PlaneID', StackName{k}];
        handles.Pipeline.(fieldname) =  ZPositions{k};
        
    end % ends looping of Groups
    handles.Current.NumberOfImageSets = length(usubRCTSRef);
end % ends first cycle calculations



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZE STACKS OF CURRENT CYCLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FileNames = cell(1, length(TextToFind));
%FileNamesWPath = cell(1, length(TextToFind));

for n=1:numGroups
    
    %%% Determines the file name of the images you want to analyze. In
    %%% contrast to the ordinary loading function, multiple columns (which
    %%% have the images of the different Z become selected
    fieldname = ['FileList', StackName{n}];
    FileList = handles.Pipeline.(fieldname);
    CurrentFileName = FileList(SetBeingAnalyzed,:);
    
    %%% obtain full path of image files
    fieldname = ['Pathname', StackName{n}];
    Pathname = handles.Pipeline.(fieldname);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% DEFINE FORMAT OF 3D INFORMATION
    %%%
    %%% Define Structure, which is at same position where in 2D
    %%% cellprofiler the Image itself would have been saved. While this
    %%% change does not affect backward compatibility with preexisting 2D
    %%% modules, it allows to reuse the interface of CP to relay
    %%% 3D information
    %%%
    %%% .Format defines Format of 3D information
    %%%
    %%% For handling 3D information it seems useful to save different 3D
    %%% information differentially. For instance direct loading of all
    %%% images could reduce the amount of spending time on loading images
    %%% at the cost of RAM. Defining only the names could allow to load
    %%% images sequentially wihtin indiviudal modules, thus reducing the
    %%% amount of RAM required for temporally storing 3D infomration.
    %%% Segmentations, might only be stored as linear indices of the full
    %%% 3D matrix to define objects of interests.
    %%%
    %%% 'NameStack' is used here to indicate:
    %%% -   path + extension were saved
    %%% -   in a cell, where
    %%%     rows indicate individual sets of stacks and
    %%%     columns individual images belonging to one set
    
    
    FileNames{n} = CurrentFileName; % for saving file names in measurments, do not add path to remain consistent w CP1
    %CurrentFileNameWPath = cellfun(@(x) fullfile(Pathname, char(x)), CurrentFileName, 'UniformOutput', 0);
    %FileNamesWPath{n}=CurrentFileNameWPath;
    handles.Pipeline.(StackName{n}).Format = 'NameStack';
    %handles.Pipeline.(StackName{n}).FileNames = FileNamesWPath{n};
    handles.Pipeline.(StackName{n}).FileNames = CurrentFileName;
    handles.Pipeline.(StackName{n}).Pathname = Pathname;
    
    fieldname = ['PlaneID', StackName{n}];
    handles.Pipeline.(StackName{n}).PlaneIDs = handles.Pipeline.(fieldname);        % add ID of initialized files
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The data with the information of the processed images will be saved
% twice in a slightly different way:
% 1)    The CP3D introdcued STACK contains full information for 3D.
% 2)    For backward compatibility IMAGE is filled with a file name of the
%       first plane of each stack.
% Both 1) and 2) follow the logic of the codde of the CP1 LoadImages module


%%%%%%% 1) + 2) COMMON REFERENCES   %%%%%%%%%%%%%%%%%

%%% First, fix feature names and the pathname
RefPathNames = cell(1,length(StackName));
RefFileNamesText = cell(1,length(StackName));
RefPathNamesText = cell(1,length(StackName));
for n = 1:length(StackName)
    RefPathNames{n} = Pathname;
    RefFileNamesText{n} = [StackName{n}];
    RefPathNamesText{n} = [StackName{n}];
end
RefFileNames = FileNames;


%%%%%%% 1)    STACK DATA    %%%%%%%%%

%%% NOTE: The structure for filenames and pathnames will be a cell array of cell arrays

%%% Since there may be several initialize modules in the pipeline which all
%%% write to the handles.Measurements.Stack.FileName field, we store
%%% filenames in an "appending" style. Here we check if any of the modules
%%% above the current module in the pipeline has written to
%%% handles.Measurements.Stack.Filenames. Then we should append the current
%%% filenames and path names to the already written ones. If this is the
%%% first module to put anything into the handles.Measurements.Stack
%%% structure, then this section is skipped and the FileNamesText fields
%%% are created with their initial entry coming from this module.

if  isfield(handles,'Measurements') && isfield(handles.Measurements,'Stack') &&...
        isfield(handles.Measurements.Stack,'FileNames') && length(handles.Measurements.Stack.FileNames) == SetBeingAnalyzed
    % Get existing file/path names. Returns a cell array of names
    StackExistingFileNamesText = handles.Measurements.Stack.FileNamesText;
    StackExistingFileNames     = handles.Measurements.Stack.FileNames{SetBeingAnalyzed};
    StackExistingPathNamesText = handles.Measurements.Stack.PathNamesText;
    StackExistingPathNames     = handles.Measurements.Stack.PathNames{SetBeingAnalyzed};
    % Append current file names to existing file names
    StackFileNamesText = cat(2,StackExistingFileNamesText,RefFileNamesText);
    StackFileNames     = cat(2,StackExistingFileNames,RefFileNames);
    StackPathNamesText = cat(2,StackExistingPathNamesText,RefPathNamesText);
    StackPathNames     = cat(2,StackExistingPathNames,RefPathNames);
else    % if no entry (eg. during first cycle), this addition became necessary to have a clean separation of the reference used by IMAGE and STOCK
    StackFileNamesText = RefFileNamesText;
    StackFileNames = RefFileNames;
    StackPathNamesText = RefPathNamesText;
    StackPathNames = RefPathNames;
end

%%% Write to the handles.Measurements.Stack structure
handles.Measurements.Stack.FileNamesText                   = StackFileNamesText;
handles.Measurements.Stack.FileNames(SetBeingAnalyzed)         = {StackFileNames};
handles.Measurements.Stack.PathNamesText                   = StackPathNamesText;
handles.Measurements.Stack.PathNames(SetBeingAnalyzed)         = {StackPathNames};



%%%%%%% 2)    IMAGE DATA    %%%%%%%%%

% obatin name of first Files
NameOfFirstFiles = cellfun(@(x) x{1}, RefFileNames, 'UniformOutput', false);

%%% Again append. Note that while this basically copies the code from
%%% above, the existing filenames in the .Image and .Stack can be
%%% different, if a loadimages module initialized into .Image thereby
%%% changeing the number of monitored Images

if  isfield(handles,'Measurements') && isfield(handles.Measurements,'Image') &&...
        isfield(handles.Measurements.Image,'FileNames') && length(handles.Measurements.Image.FileNames) == SetBeingAnalyzed
    % Get existing file/path names. Returns a cell array of names
    ImageExistingFileNamesText = handles.Measurements.Image.FileNamesText;
    ImageExistingFileNames     = handles.Measurements.Image.FileNames{SetBeingAnalyzed};
    ImageExistingPathNamesText = handles.Measurements.Image.PathNamesText;
    ImageExistingPathNames     = handles.Measurements.Image.PathNames{SetBeingAnalyzed};
    % Append current file names to existing file names
    ImageFileNamesText = cat(2,ImageExistingFileNamesText,RefFileNamesText);
    ImageFileNames     = cat(2,ImageExistingFileNames,NameOfFirstFiles);
    ImagePathNamesText = cat(2,ImageExistingPathNamesText,RefPathNamesText);
    ImagePathNames     = cat(2,ImageExistingPathNames,RefPathNames);
else    % if no entry (eg. during first cycle), this addition became necessary to have a clean separation of the reference used by IMAGE and STOCK
    ImageFileNamesText = RefFileNamesText;
    ImageFileNames = NameOfFirstFiles;
    ImagePathNamesText = RefPathNamesText;
    ImagePathNames = RefPathNames;
end


%%% Write to the handles.Measurements.Image structure
handles.Measurements.Image.FileNamesText                   = ImageFileNamesText;
handles.Measurements.Image.FileNames(SetBeingAnalyzed)         = {ImageFileNames};
handles.Measurements.Image.PathNamesText                   = ImagePathNamesText;
handles.Measurements.Image.PathNames(SetBeingAnalyzed)         = {ImagePathNames};

end