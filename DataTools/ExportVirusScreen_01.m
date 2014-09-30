function ExportVirusScreen_01(handles)
global bsHack

% Help for the ExportVirusScreen tool:
% Category: Data Tools
%
% SHORT DESCRIPTION:
% Exports center locations of objects. Specialty function for creating a
% locations list for microscopy image acquisition of gridded spots.
% *************************************************************************
% Useful for creating a locations list for microscope.

%%% Ask the user to choose the file from which to extract
%%% measurements. The window opens in the default output directory.
%%% [BS] But only if the current handles does not contain any measurements.
%%% (for instance, when called directly from the VirusScreen module).

if ~isfield(handles,'Measurements')

    [RawFileName, RawPathname] = uigetfile(fullfile(handles.Current.DefaultOutputDirectory,'.','*.mat'),'Select the raw measurements file');
    
    %%% Allows canceling.
    if RawFileName == 0
        return
    end
    load(fullfile(RawPathname, RawFileName));

    if ~exist('handles','var')
        CPerrordlg('This is not a CellProfiler output file.');
        return
    end

    %%% Quick check if it seems to be a CellProfiler file or not
    if ~isfield(handles,'Measurements')
        CPerrordlg('The selected file does not contain any measurements.')
        return
    end

    %%% Extract the fieldnames of measurements from the handles structure.
    MeasFieldnames = fieldnames(handles.Measurements);
    for i=1:length(handles.Measurements)
        if strcmp(MeasFieldnames{i},'Image')
            MeasFieldnames(i) = [];
        end
        if isempty(MeasFieldnames)
            CPerrordlg('The output file you have chosen does not have location measurements.');
            return
        end
    end
    [Selection, ok] = listdlg('ListString',MeasFieldnames,'ListSize', [300 200],...
        'Name','Select measurement',...
        'PromptString','Which object do you want to export locations from?',...
        'CancelString','Cancel',...
        'SelectionMode','single');

    if ok == 0
        return
    end

    ObjectTypename = MeasFieldnames{Selection};    
    
    if ~isfield(handles.Measurements.(ObjectTypename),'VirusInfection')
        CPerrordlg('The object you have chosen does not have any VirusScreen data.');
        return
    end
    
else

    %loop through all measurement objects searching for VirusInfection
    %data.
    MeasFieldnames = fieldnames(handles.Measurements);
    %CPerrordlg(['[BS] Objects in Measurements: ',MeasFieldnames{:}, ' - ',num2str(length(MeasFieldnames))]);
    for ix=1:length(MeasFieldnames)
        if ~strcmp(MeasFieldnames{ix},'Image')
            %CPerrordlg(['[BS] Found an object that is not Image!']);
            ObjectFieldnames = fieldnames(handles.Measurements.(MeasFieldnames{ix}));
            for ii=1:length(ObjectFieldnames)
                if strcmp(ObjectFieldnames{ii},'VirusInfection')
                    ObjectTypename = MeasFieldnames{ix};
                    %CPhelpdlg(['[BS] Exported the data from your current ''',ObjectTypename,''' object.']);
                end
            end
        end
    end
    
    
    % If no VirusInfection data is found, repeat request for Data file...
    % (ugly copypasting of code...)
    if ~exist('ObjectTypename')
        %CPerrordlg('[BS] Your current data does not contain any VirusScreen data. Please select an old data file.');

        [RawFileName, RawPathname] = uigetfile(fullfile(handles.Current.DefaultOutputDirectory,'.','*.mat'),'Select the raw VirusScreen output file')

        %%% Allows canceling.
        if RawFileName == 0
            return
        end
        load(fullfile(RawPathname, RawFileName));

        if ~exist('handles','var')
            CPerrordlg('This is not a CellProfiler output file.');
            return
        end

        %%% Quick check if it seems to be a CellProfiler file or not
        if ~isfield(handles,'Measurements')
            CPerrordlg('The selected file does not contain any measurements.')
            return
        end

        %%% Extract the fieldnames of measurements from the handles structure.
        MeasFieldnames = fieldnames(handles.Measurements);
        for i=1:length(handles.Measurements)
            if strcmp(MeasFieldnames{i},'Image')
                MeasFieldnames(i) = [];
            end
            if isempty(MeasFieldnames)
                CPerrordlg('The output file you have chosen does not have location measurements.');
                return
            end
        end
        [Selection, ok] = listdlg('ListString',MeasFieldnames,'ListSize', [300 200],...
            'Name','Select measurement',...
            'PromptString','Which object do you want to export locations from?',...
            'CancelString','Cancel',...
            'SelectionMode','single');

        if ok == 0
            return
        end

        ObjectTypename = MeasFieldnames{Selection};    

        if ~isfield(handles.Measurements.(ObjectTypename),'VirusInfection')
            CPerrordlg('The object you have chosen does not have any VirusScreen data.');
            return
        end
    end
    
    
end

format short;


intVirusScreenModuleNumber = 0;
intDoPostAnlysisNumber = 0;
strDoPostAnalysis = '';
intDiscardImageFeatureNumber = 0;
strListOfDiscardedImages = '';
iiold = 0;
ii = 0;

if isfield(handles.Settings, 'VariableNames')
    intVirusScreenModuleNumber = strmatch('VirusScreen', handles.Settings.ModuleNames);
    
% [BS 2006-06-17] Standardize postanalysis again

    strDoPostAnalysis = 'Yes';

%     intDoPostAnlysisNumber = strmatch('strDoPostAnalysis', handles.Settings.VariableNames{intVirusScreenModuleNumber});
%     if not(intDoPostAnlysisNumber == 0)
%         strDoPostAnalysis = handles.Settings.VariableValues{intVirusScreenModuleNumber,intDoPostAnlysisNumber};
%     else
%         strDoPostAnalysis = 'No';
%     end
    intDiscardImageFeatureNumber = strmatch('DiscardImage', handles.Measurements.(ObjectTypename).VirusInfectionFeatures);
end

%[BS] Raw output file
if length(handles.Measurements.Image.FileNames{1}{1}) > 13
    strFileSuffix = handles.Measurements.Image.FileNames{1}{1}(1:end-13);
else
    strFileSuffix = ObjectTypename;
end

tempding = strfind(strFileSuffix, '/');
if not(isempty(tempding))
    strFileSuffix = strFileSuffix(tempding(end)+1:end);
end
filename = [Datestr(now,29),'_',strrep(Datestr(now,13),':',''),'_',strFileSuffix,'_VirusScreen_raw.csv'];    
fid = fopen(fullfile(handles.Current.DefaultOutputDirectory,filename),'w');

if fid == -1
    fid = fopen(fullfile(RawPathname,filename),'w');
    CPmsgbox(sprintf('Could not load from here: %s%s',handles.Current.DefaultOutputDirectory,filename));
    if fid == -1
        CPerrordlg(sprintf('Could also not load from here: %s%s',RawPathname,filename));        
        %CPerrordlg(sprintf('Cannot create the output file %s. There might be another program using a file with the same name.',filename));
        return
    end
end

fprintf(fid,'FileName;');
for i = 1:length(handles.Measurements.(ObjectTypename).VirusInfectionFeatures)
    fprintf(fid,'%s;',handles.Measurements.(ObjectTypename).VirusInfectionFeatures{i});
end
fprintf(fid,'\n');

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

for i = 1:length(handles.Measurements.(ObjectTypename).VirusInfection)
    fprintf(fid,'%s;',handles.Measurements.Image.FileNames{i}{1});
    for iix = 1:length(handles.Measurements.(ObjectTypename).VirusInfection{i}(1,:))
        fprintf(fid,'%1.3f;',handles.Measurements.(ObjectTypename).VirusInfection{i}(1,iix));
    end
    
    if strcmp(strDoPostAnalysis, 'Yes') || length(strmatch(ModuleName, 'VirusScreen_Auto')) > 0
        if handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intDiscardImageFeatureNumber) == 1
            strListOfDiscardedImages = strcat(strListOfDiscardedImages, ' - C',num2str(i+1));
        end
        strFileName = handles.Measurements.Image.FileNames{i}{1};
        strWellID = strFileName(end-11:end-9);
        % write the sum() equation at every last well image
        % either, next image is new well, or at last image :)
        % include the collected discarded images string.
        if i == length(handles.Measurements.(ObjectTypename).VirusInfection) || length(findstr(handles.Measurements.Image.FileNames{i}{1}, handles.Measurements.Image.FileNames{i+1}{1}(1:end-9))) == 0 
            iiold = ii;
            ii = i;
            fprintf(fid,'%s;', strcat('=SUM(B', num2str(iiold+2), ':B', num2str(ii+1), ')', strListOfDiscardedImages));
            fprintf(fid,'%s;', strcat('=SUM(C', num2str(iiold+2), ':C', num2str(ii+1), ')', strrep(strListOfDiscardedImages, 'B', 'C')));
            strListOfDiscardedImages = '';
        end
    end
    
    fprintf(fid,'\n');
end    

fclose(fid);

if not(bsHack)
    CPhelpdlg(['[BS] Saved ',handles.Current.DefaultOutputDirectory,'\',filename]);
end