function ExportVirusScreenLocalDensity_01(handles)
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

        if ~isfield(handles.Measurements.(ObjectTypename),'VirusScreenLDCorrectedPlateData')
            CPerrordlg('The object you have chosen does not have any Local Density VirusScreen data.');
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
filename = [datestr(now,29),'_',strrep(datestr(now,13),':',''),'_',strFileSuffix,'_localdensity.csv'];    
fid = fopen(fullfile(handles.Current.DefaultOutputDirectory,filename),'w');

try
    if fid == -1
        fid = fopen(fullfile(RawPathname,filename),'w');
        CPmsgbox(sprintf('Could not load from here: %s%s',handles.Current.DefaultOutputDirectory,filename));
        if fid == -1
            CPerrordlg(sprintf('Could also not load from here: %s%s',RawPathname,filename));        
            %CPerrordlg(sprintf('Cannot create the output file %s. There might be another program using a file with the same name.',filename));
            return
        end
    end
catch
    %[BS] try on windows PC... hardcode HACK
    fid = fopen(fullfile('C:\Documents and Settings\imsb\Desktop',filename),'w');
end

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

cellLDCorrDataFeatures = handles.Measurements.(ObjectTypename).VirusScreenLDCorrectedPlateDataFeatures;

for i = 1:length(cellLDCorrDataFeatures)
    fprintf(fid,'%s;',cellLDCorrDataFeatures{i});
end
fprintf(fid,'\n');

cellLDCorrData = handles.Measurements.(ObjectTypename).VirusScreenLDCorrectedPlateData;

for ii = 1:size(handles.Measurements.Nuclei.VirusScreenLDCorrectedPlateData,1)
    fprintf(fid,'%s;',cellLDCorrData{ii,1});
    fprintf(fid,'%s;',cellLDCorrData{ii,2});
    fprintf(fid,'%s;',cellLDCorrData{ii,3});    
    for iix = 4:size(cellLDCorrData,2)
        fprintf(fid,'%1.3f;',cellLDCorrData{ii,iix});
    end
    fprintf(fid,'\n');
end    

fclose(fid);

if not(bsHack)
    CPhelpdlg(['[BS] Saved Local Density data to ',handles.Current.DefaultOutputDirectory,'\',filename]);
end