function handles = MinimalInitialize3DStack(handles)

% Help for the Load Images module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% This is an experimental module for testing a 3D workflow with CP
% 
% TO DO:
% x Add security checks ensuring that images from ordinary 2d Load..Module
%   and Initialize3D refer to the same image set
% x See how the create Batch files would work with this module
% x Unlock multiple inputs
% x support for defining whether input is exact match or regexp
% x Internally check that image sets from different stack-specifiers refer 
%   to the same position
% x check input variables
% x add support for loading all Z-stacks, eg. input 'all' instead of (1:15)
% x add offset and step size
% x fix display

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

%textVAR03 = Which Z positions do you want to load?
%defaultVAR03 = 1:10
ZPositions{1} = eval(char(handles.Settings.VariableValues{CurrentModuleNum,3}));

%pathnametextVAR04 = Enter the path name to the folder where the images to be loaded are located. Type period (.) for default image folder.
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,4});




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
for n = 1:1
    if ~strcmp(TextToFind{n}, '/') && ~strcmp(StackName{n}, '/')
        tmp1{end+1} = TextToFind{n};
        tmp2{end+1} = StackName{n};
    end
end
TextToFind = tmp1;
StackName = tmp2;
clear tmp1; clear tmp2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if SetBeingAnalyzed == 1
    %%% Specify Directory of Files
    %%% Get the pathname and check that it exists
    if strncmp(Pathname,'.',1)
        if length(Pathname) == 1
            Pathname = handles.Current.DefaultImageDirectory;
        else
            Pathname = fullfile(handles.Current.DefaultImageDirectory,Pathname(2:end));
        end
    end
    
    SpecifiedPathname = Pathname;
    if ~exist(SpecifiedPathname,'dir')
        error(['Image processing was canceled in the ', ModuleName, ' module because the directory "',SpecifiedPathname,'" does not exist. Be sure that no spaces or unusual characters exist in your typed entry and that the pathname of the directory begins with /.'])
    end
    
    numUniqueElements = NaN(1,length(TextToFind));
    
    % For each of the specified elements, extract Filenames of one stack
    for k=1: length(TextToFind)
        
        FileList = CPretrievemediafilenames(SpecifiedPathname,char(TextToFind{k}),'', 'Regular','Image')';
        %%% Checks whether any files are left.
        if isempty(FileList)
            error(['Image processing was canceled in the ', ModuleName, ' module because there are no image files with the text "', TextToFind{k}, '" in the chosen directory (or subdirectories, if you requested them to be analyzed as well).'])
        end
        
        % extract metadata of all files fitting the description of the stacks
        [intRow, intColumn, ~, intTimepoint, intZPos, intAction, intSite] = ...
            cellfun(@filterimagenamedataWithZ,FileList ,'UniformOutput',false);
        
        % Now define, which images belong together: Group by combination of
        % Row, Column, Timepoint, Site and Action
        intRCTSA= [cell2mat(intRow) cell2mat(intColumn) cell2mat(intTimepoint) cell2mat(intSite) cell2mat(intAction)];
        maxintRCTSA = max(intRCTSA,[],1); % Find maxima to create subindex
        subRCTSA = sub2ind2(maxintRCTSA,intRCTSA);
        intZPos = cell2mat(intZPos);
        
        % Filter dataset for first Z position specified, this should yield a
        % unique list of combinations of Rows, Columns, Timepoint, Site and
        % Action (which will be checked later)
        f = intZPos == ZPositions{k}(1,1);
        fsubRCTSA = subRCTSA(f,1);
        
        % Monitor Amount of unique elements for each different Category of
        % stack (e.g.: StackBlue, StackGreen)..
        
        numUniqueElements(1,k) = length(unique(fsubRCTSA));
        
        % Check if identification of images is unambiguous
        if  numUniqueElements(1,k)~= length(fsubRCTSA)
            error(['Image processing was canceled in the ', ModuleName, ' module because there are "', TextToFind{k}, '" is not unambiguous for a given Rows, Columns, Site, Timepoint combination.'])
        end
        
        % If number of Stacks for each category of stacks is equal,
        % use this as the number of image sets, if not the specified inputs
        % refer to a different number of
        
        handles.Current.NumberOfImageSets = numUniqueElements(1,k);
        
        % Read file names belonging to a given set of stacks. Rows will
        % indicate different Stacks, whereas Columns indicate different images
        % belonging to one Stack
        dimImageStacks = size(fsubRCTSA,1);
        dimImagesPerStack = size(ZPositions{k},2);
        
        StackFileList = cell(dimImageStacks,dimImagesPerStack);
        for l=1:dimImageStacks
            for r = 1: dimImagesPerStack
                currentfsubRCTSA = fsubRCTSA(l,1);
                f = (subRCTSA == currentfsubRCTSA) & (intZPos==ZPositions{k}(1,r));
                StackFileList{l,r}=FileList(f,1);
            end
        end
        
        % Output to handles Structure, analogously to conventional loading
        % function
        fieldname = ['FileList', StackName{k}];
        handles.Pipeline.(fieldname) = StackFileList;
        fieldname = ['Pathname', StackName{k}];
        handles.Pipeline.(fieldname) = SpecifiedPathname;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZE STACKS OF CURRENT CYCLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FileNames = cell(1, length(TextToFind));


for n=1:length(TextToFind)
    
    %%% Determines the file name of the images you want to analyze. In
    %%% contrast to the ordinary loading function, multiple columns (which
    %%% have the images of the different Z become selected
    fieldname = ['FileList', StackName{n}];
    FileList = handles.Pipeline.(fieldname);
    CurrentFileName = FileList(SetBeingAnalyzed,:);
    
    %%% obtain full path of image files
    fieldname = ['Pathname', StackName{n}];
    Pathname = handles.Pipeline.(fieldname);
    
    CurrentFileName = cellfun(@(x) fullfile(Pathname, char(x)), CurrentFileName, 'UniformOutput', 0);
   
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
    %%% -   full filenames (path + extension) were saved
    %%% -   in a cell, where
    %%%     rows indicate individual sets of stacks and
    %%%     columns individual images belonging to one set
    
    FileNames{n}=CurrentFileName;
    handles.Pipeline.(StackName{n}).FileNames = FileNames{n};
    handles.Pipeline.(StackName{n}).Format = 'NameStack';
end



%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure('','NarrowText',ThisModuleFigureNumber)
    end
    for n = 1:length(StackName)
        %%% Activates the appropriate figure window.
        currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
        if iscell(StackName)
            TextString = [StackName{n},': ',handles.Pipeline.(StackName{n}).FileNames];
        else
            TextString = [StackName,': ',handles.Pipeline.(StackName).FileNames];
        end
        %TO DO: fix so that all the files are shown there are some, which
        %only become visible after manually resizing the window
        
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','left','string',TextString,'position',[.05 0.85-(n-1)*.15 .95 .1],'BackgroundColor',[.7 .7 .9])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTE: The structure for filenames and pathnames will be a cell array of cell arrays

%%% First, fix feature names and the pathname
PathNames = cell(1,length(StackName));
FileNamesText = cell(1,length(StackName));
PathNamesText = cell(1,length(StackName));
for n = 1:length(StackName)
    PathNames{n} = Pathname;
    FileNamesText{n} = [StackName{n}];
    PathNamesText{n} = [StackName{n}];
end

%%% Since there may be several load/initialize/save modules in the pipeline which all
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
    ExistingFileNamesText = handles.Measurements.Stack.FileNamesText;
    ExistingFileNames     = handles.Measurements.Stack.FileNames{SetBeingAnalyzed};
    ExistingPathNamesText = handles.Measurements.Stack.PathNamesText;
    ExistingPathNames     = handles.Measurements.Stack.PathNames{SetBeingAnalyzed};
    % Append current file names to existing file names
    FileNamesText = cat(2,ExistingFileNamesText,FileNamesText);
    FileNames     = cat(2,ExistingFileNames,FileNames);
    PathNamesText = cat(2,ExistingPathNamesText,PathNamesText);
    PathNames     = cat(2,ExistingPathNames,PathNames);
end

%%% Write to the handles.Measurements.Stack structure
handles.Measurements.Stack.FileNamesText                   = FileNamesText;
handles.Measurements.Stack.FileNames(SetBeingAnalyzed)         = {FileNames};
handles.Measurements.Stack.PathNamesText                   = PathNamesText;
handles.Measurements.Stack.PathNames(SetBeingAnalyzed)         = {PathNames};

end