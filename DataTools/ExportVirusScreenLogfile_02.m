function ExportVirusScreenLogfile_02(handles)
% Help for the ExportVirusScreenLogfile_02 tool:
% Category: Data Tools
%
% SHORT DESCRIPTION:
% Exports center locations of objects. Specialty function for creating a
% locations list for microscopy image acquisition of gridded spots.
% *************************************************************************
% Useful for creating a locations list for microscope.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Susan Ma
%   Wyman Li
%
% Website: http://www.cellprofiler.org
%
% $Revision: 2950 $

%%% Ask the user to choose the file from which to extract
%%% measurements. The window opens in the default output directory.
%%% [BS] But only if the current handles does not contain any measurements.
%%% (for instance, when called directly from the VirusScreen module).

RawPathname = [];

if ~isfield(handles,'Measurements')

    [RawFileName, RawPathname] = uigetfile(fullfile(handles.Current.DefaultOutputDirectory,'.','*.mat'),'Select the raw measurements file');
    
    %%% Allows canceling.
    if RawFileName == 0
        return
    end
    load(fullfile(RawPathname, RawFileName));

    if ~exist('handles','var')
        CPwarndlg('This is not a CellProfiler output file.');
        return
    end

    %%% Quick check if it seems to be a CellProfiler file or not
    if ~isfield(handles,'Measurements')
        CPwarndlg('The selected file does not contain any measurements.')
        return
    end

    %%% Extract the fieldnames of measurements from the handles structure.
    MeasFieldnames = fieldnames(handles.Measurements);
    for i=1:length(handles.Measurements)
        if strcmp(MeasFieldnames{i},'Image')
            MeasFieldnames(i) = [];
        end
        if isempty(MeasFieldnames)
            CPwarndlg('The output file you have chosen does not have location measurements.');
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
        CPwarndlg('The object you have chosen does not have any VirusScreen data.');
        return
    end
    
else

    %loop through all measurement objects searching for VirusInfection
    %data.
    MeasFieldnames = fieldnames(handles.Measurements);
    %CPwarndlg(['[BS] Objects in Measurements: ',MeasFieldnames{:}, ' - ',num2str(length(MeasFieldnames))]);
    for ix=1:length(MeasFieldnames)
        if ~strcmp(MeasFieldnames{ix},'Image')
            %CPwarndlg(['[BS] Found an object that is not Image!']);
            ObjectFieldnames = fieldnames(handles.Measurements.(MeasFieldnames{ix}));
            for ii=1:length(ObjectFieldnames)
                if strcmp(ObjectFieldnames{ii},'VirusInfection')
                    ObjectTypename = MeasFieldnames{ix};
                    %CPhelpdlg(['[BS] Exported the logfile from your current ''',ObjectTypename,''' object.']);
                end
            end
        end
    end
    
    
    % If no VirusInfection data is found, repeat request for Data file...
    % (ugly copypasting of code...)
    if ~exist('ObjectTypename')
        %CPwarndlg('[BS] Your current data does not contain any VirusScreen data. Please select an old data file.');

        [RawFileName, RawPathname] = uigetfile(fullfile(handles.Current.DefaultOutputDirectory,'.','*.mat'),'Select the raw measurements file');

        %%% Allows canceling.
        if RawFileName == 0
            return
        end
        load(fullfile(RawPathname, RawFileName));

        if ~exist('handles','var')
            CPwarndlg('This is not a CellProfiler output file.');
            return
        end

        %%% Quick check if it seems to be a CellProfiler file or not
        if ~isfield(handles,'Measurements')
            CPwarndlg('The selected file does not contain any measurements.')
            return
        end

        %%% Extract the fieldnames of measurements from the handles structure.
        MeasFieldnames = fieldnames(handles.Measurements);
        for i=1:length(handles.Measurements)
            if strcmp(MeasFieldnames{i},'Image')
                MeasFieldnames(i) = [];
            end
            if isempty(MeasFieldnames)
                CPwarndlg('The output file you have chosen does not have location measurements.');
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
            CPwarndlg('The object you have chosen does not have any VirusScreen data.');
            return
        end
    end
    
    
end

if length(handles.Measurements.Image.FileNames{1}{1}) > 13
    strFileSuffix = handles.Measurements.Image.FileNames{1}{1}(1:end-13);
else
    strFileSuffix = ObjectTypename;
end
tempding = strfind(strFileSuffix, '/');
if not(isempty(tempding))
    strFileSuffix = strFileSuffix(tempding(end)+1:end);
end

filename = [Datestr(now,29),'_',strrep(Datestr(now,13),':',''),'_',strFileSuffix,'_logfile.txt'];

% save logfile next to .mat outputfile if this has been hand picked.
if not(isempty(RawPathname))
    fid = fopen(fullfile(RawPathname,filename),'w');
    disp(fullfile(RawPathname,filename))
else
% otherwise save in current default output directory    
    fid = fopen(fullfile(handles.Current.DefaultOutputDirectory,filename),'w');
    disp(fullfile(handles.Current.DefaultOutputDirectory,filename))
end

if fid == -1
    CPwarndlg(sprintf('Cannot create the output file %s. There might be another program using a file with the same name.',fullfile(handles.Current.DefaultOutputDirectory,filename)));
    return
end

intVirusScreenModuleNumber = 0;
intDiscardImageFeatureNumber = 0;
intOutOfFocusFeatureNumber = 0;
intDiscardImageFeatureNumber = 0;
intTotalImages = 0;
intSettingsLength = 0;
intOutOfFocusCutoff = inf;
intDiscarded = 0;
intDiscardedControlImages = 0;
intTotalControlImages = 0;
intTwoSingleGaussians = 0;
intOneSingleGaussians = 0;
intOutOfFocusImages = 0;
intGaussianAndTooLittleNuclei = 0;
intIdentifyPrimAutomaticModuleNumber=0;

if isfield(handles.Measurements.(ObjectTypename),'VirusScreenDiscardedControls')
    intDiscardedControlImages = handles.Measurements.(ObjectTypename).VirusScreenDiscardedControls{1}(1,1);
end
if isfield(handles.Measurements.(ObjectTypename),'VirusScreenTotalControls')
    intTotalControlImages = handles.Measurements.(ObjectTypename).VirusScreenTotalControls{1}(1,1);
end

strLDCtrlWells = [];
cellLDCorrPlateData = {};
intCtrlII = [];
intCtrlCN = [];

if isfield(handles.Measurements.(ObjectTypename),'VirusScreenLDCorrectedPlateData')
    if not(isempty(handles.Measurements.(ObjectTypename).VirusScreenLDCorrectedPlateData))
        intLDPlateDataCtrlIInumbers = strmatch('CtrlII', handles.Measurements.(ObjectTypename).VirusScreenLDCorrectedPlateDataFeatures);
        intLDPlateDataCtrlCNnumbers = strmatch('CtrlCN', handles.Measurements.(ObjectTypename).VirusScreenLDCorrectedPlateDataFeatures);
        cellLDCorrPlateData = handles.Measurements.(ObjectTypename).VirusScreenLDCorrectedPlateData;
        intCtrlII = [cellLDCorrPlateData{:,intLDPlateDataCtrlIInumbers}];
        intCtrlCN = [cellLDCorrPlateData{:,intLDPlateDataCtrlCNnumbers}];
    end
end

try
    intVirusScreenModuleNumber = strmatch('VirusScreen', handles.Settings.ModuleNames);
    intLoadImagesModuleNumber = strmatch('LoadImages', handles.Settings.ModuleNames);
    intRescaleIntensityModuleNumber = strmatch('RescaleIntensity', handles.Settings.ModuleNames);
    intIdentifyPrimAutomaticModuleNumber = strmatch('IdentifyPrimAutomatic', handles.Settings.ModuleNames);    
    intDiscardImageFeatureNumber = strmatch('DiscardImage', handles.Measurements.(ObjectTypename).VirusInfectionFeatures);
    intFirstSingleFeatureNumber = strmatch('FirstSingleGaussian', handles.Measurements.(ObjectTypename).VirusInfectionFeatures);
    intSecondSingleFeatureNumber = strmatch('SecondSingleGaussian', handles.Measurements.(ObjectTypename).VirusInfectionFeatures);
    intOutOfFocusFeatureNumber = strmatch('OutOfFocus', handles.Measurements.(ObjectTypename).VirusInfectionFeatures);
    intTotalImages = length(handles.Measurements.Image.FileNames);
    
    intOutOfFocusVariableNumber = strmatch('intOutOfFocusCutoff',handles.Settings.VariableNames{intVirusScreenModuleNumber});
    if not(isempty(intOutOfFocusVariableNumber))
        intOutOfFocusCutoff = str2double(handles.Settings.VariableValues{intVirusScreenModuleNumber,intOutOfFocusVariableNumber});
    end
    if isfield(handles.Settings, 'VariableNames')
        intSettingsLength = double(length(handles.Settings.VariableNames{intVirusScreenModuleNumber}));
    end
    
    intLDCtrlWellsVariableNumber = strmatch('strLocalDensityCtrlWells',handles.Settings.VariableNames{intVirusScreenModuleNumber});
    if not(isempty(intLDCtrlWellsVariableNumber))
        strLDCtrlWells = char(handles.Settings.VariableValues{intVirusScreenModuleNumber,intLDCtrlWellsVariableNumber});
    end

end

fprintf(fid, 'BS CELLPROFILER VIRUSSCREEN LOGFILE\n\n');
fprintf(fid, 'Run Started on: \t%s\n', handles.Current.TimeStarted);
fprintf(fid, 'Log Created on: \t%s\n', Datestr(now,0));
if exist('RawFileName')
    fprintf(fid, 'Raw data file: \t%s\n', [RawPathname, RawFileName]);
else
    fprintf(fid, 'Output folder: \t%s\n', handles.Current.DefaultOutputDirectory);    
end
fprintf(fid, 'Image folder: \t%s\n', handles.Current.DefaultImageDirectory);

fprintf(fid, '\nPOST ANALYSIS OVERVIEW\n');
for i = 1:length(handles.Measurements.(ObjectTypename).VirusInfection)
    if handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intDiscardImageFeatureNumber) == 1
        intDiscarded = intDiscarded + 1;
        if handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intOutOfFocusFeatureNumber) > intOutOfFocusCutoff
            intOutOfFocusImages = intOutOfFocusImages + 1;
        else
            intGaussianAndTooLittleNuclei = intGaussianAndTooLittleNuclei + 1;
        end
    end
    if not(isempty(intFirstSingleFeatureNumber)) && (handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intFirstSingleFeatureNumber) == 1 || handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intSecondSingleFeatureNumber) == 1) 
        if (handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intFirstSingleFeatureNumber) == 1 && handles.Measurements.(ObjectTypename).VirusInfection{i}(1,intSecondSingleFeatureNumber) == 1)         
            intTwoSingleGaussians = intTwoSingleGaussians + 1;
        else
            intOneSingleGaussians = intOneSingleGaussians + 1;
        end
    end    
end

fprintf(fid, 'Total images analyzed     \t = %d\n',intTotalImages);
fprintf(fid, '  Total discrd imgs       \t = %d (%1.1f%%)\n',intDiscarded, single((intDiscarded/intTotalImages)*100));
fprintf(fid, '  Total ctrl imgs         \t = %d\n',intTotalControlImages);
fprintf(fid, '  Total discrd ctrl imgs  \t = %d\n\n',intDiscardedControlImages);

fprintf(fid, '  Local density ctrl wells\t = %s\n',strLDCtrlWells);
fprintf(fid, '  Mean ctrl well total    \t = %.0f\n',round(intCtrlCN));
fprintf(fid, '  Mean ctrl well ii       \t = %1.3f\n',intCtrlII);

fprintf(fid, '\nUSED SETTINGS\n');
fprintf(fid, 'Module #%.0f: %s\n',intVirusScreenModuleNumber,handles.Settings.ModuleNames{intVirusScreenModuleNumber});
for i = 1:intSettingsLength
    fprintf(fid, '%s = %s\n', [char(handles.Settings.VariableNames{intVirusScreenModuleNumber}(1,i)), blanks(30 - length(char(handles.Settings.VariableNames{intVirusScreenModuleNumber}(1,i))))], char(handles.Settings.VariableValues{intVirusScreenModuleNumber,i}));
end

fprintf(fid, '\nModule #%.0f: %s\n',intLoadImagesModuleNumber,handles.Settings.ModuleNames{intLoadImagesModuleNumber});
fprintf(fid, '%s = %s\n', ['Mode', blanks(30 - length('Mode'))], char(handles.Settings.VariableValues{intLoadImagesModuleNumber,1}));
fprintf(fid, '%s = %s\n', [char(handles.Settings.VariableValues{intLoadImagesModuleNumber,3}), blanks(30 - length(char(handles.Settings.VariableValues{intLoadImagesModuleNumber,3})))], char(handles.Settings.VariableValues{intLoadImagesModuleNumber,2}));
fprintf(fid, '%s = %s\n', [char(handles.Settings.VariableValues{intLoadImagesModuleNumber,5}), blanks(30 - length(char(handles.Settings.VariableValues{intLoadImagesModuleNumber,5})))], char(handles.Settings.VariableValues{intLoadImagesModuleNumber,4}));
if not(isempty(char(handles.Settings.VariableValues{intLoadImagesModuleNumber,15})))
    fprintf(fid, '%s = %s\n', ['Control', blanks(30 - length('Control'))], char(handles.Settings.VariableValues{intLoadImagesModuleNumber,15}));
end

fprintf(fid, '\nModule #%.0f: %s\n',intRescaleIntensityModuleNumber,handles.Settings.ModuleNames{intRescaleIntensityModuleNumber});
fprintf(fid, '%s = %s\n', ['New max intensity', blanks(30 - length('New max intensity'))], char(handles.Settings.VariableValues{intRescaleIntensityModuleNumber,5}));
try
    ImageName = 'OrigGreen';%char(handles.Settings.VariableValues{intVirusScreenModuleNumber,2});
    fieldname = ['MinPixelValue', ImageName];
    fprintf(fid, '%s = %1.6f\n', [fieldname, blanks(30 - length(fieldname))], double(handles.Pipeline.(fieldname)));
    fieldname = ['MaxPixelValue', ImageName];
    fprintf(fid, '%s = %1.6f\n', [fieldname, blanks(30 - length(fieldname))], double(handles.Pipeline.(fieldname)));
    fieldname = ['OrigMaxPixelValue', ImageName];
    fprintf(fid, '%s = %1.6f\n', [fieldname, blanks(30 - length(fieldname))], double(handles.Pipeline.(fieldname)));
end

fprintf(fid, '\nModule #%.0f: %s\n',intIdentifyPrimAutomaticModuleNumber,handles.Settings.ModuleNames{intIdentifyPrimAutomaticModuleNumber});
fprintf(fid, 'Nucleus diameter range    \t = [%s]\n', char(handles.Settings.VariableValues{intIdentifyPrimAutomaticModuleNumber,3}));
fprintf(fid, 'Method                    \t = %s\n', char(handles.Settings.VariableValues{intIdentifyPrimAutomaticModuleNumber,7}));
fprintf(fid, 'Threshold range           \t = [%s]\n', char(handles.Settings.VariableValues{intIdentifyPrimAutomaticModuleNumber,9}));

fclose(fid);