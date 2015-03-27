function handles = VirusScreen_Cluster_03(handles)

warning off all

% Help for the VirusScreen_Cluster_02 module:
% Category: Other
%
% SHORT DESCRIPTION:
%
% DOES ALL GAUSSIAN FITTING (SINGLE AND DOUBLE) ON ALL INTENSITY
% MEASUREMENTS, INCLUDING LOG10 TRANSFORMATION!!!! (VERSION '_3')
% *************************************************************************
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum] = CPwhichmodule(handles);

%textVAR01 = Which object would you like to use for the infection measurement?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Which image contains the virus signal measurements?
%infotypeVAR02 = imagegroup
%inputtypeVAR02 = popupmenu
strImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which regular expression matches your no-virus controls? (examples are: "_H\d\d" for row H, or "_[A-Z]01" for column 1)?
%defaultVAR03 = _[A-Z]01
strRegexpOfNoVirusCtrls = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Are your no-virus control images infected?
%choiceVAR04 = No
%choiceVAR04 = Yes
%inputtypeVAR04 = popupmenu
strInfectedControls = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%%VariableRevisionNumber = 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PREPARATION OF VARIABLES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    intSetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

    if (intSetBeingAnalyzed == 1)
        handles.Settings.VariableNames{CurrentModuleNum} = {
            'strObjectName',...
            'strImageName',...
            'strRegexpOfNoVirusCtrls',...
            'strInfectedControls',...
        };
    end

    strIntensityMeasurementList = '1,2,3,4,5,6,7,8,9,10,11';
    matIntensityMeasurements = str2num(strrep(strIntensityMeasurementList, ',', ' '));
    
    intSDsForThreshold = 3;

    strIntensityFieldName = ['Intensity_', strImageName];
    intTotalCells = handles.Measurements.Image.ObjectCount{intSetBeingAnalyzed}(1,1);

    IntensityMeasurements = handles.Measurements.(strObjectName).(strIntensityFieldName){intSetBeingAnalyzed}; %(:,intFirstIntensityMeasurement);

    matIntensityThresholds = [];
    cellGaussians = cell(2,length(matIntensityMeasurements), 3);

    if ~isfield(handles.Measurements.(strObjectName),'VirusScreenGaussians')
        handles.Measurements.(strObjectName).VirusScreenGaussians = {};
    end
    if ~isfield(handles.Measurements.(strObjectName),'VirusScreenThresholds')
        handles.Measurements.(strObjectName).VirusScreenThresholds = {};
    end

    try
        if intTotalCells > 1        
            % Do all one and two gaussian fittings
            for iNumberOfGaussians = 1:2
                for iIntensity = 1:11%length(matIntensityMeasurements)
                    matCurrentIntensityMeasurements = IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity));
                    matCurrentIntensityMeasurements = log10(matCurrentIntensityMeasurements);
                    matCurrentIntensityMeasurements(isinf(matCurrentIntensityMeasurements))=[];
                    matCurrentIntensityMeasurements(isnan(matCurrentIntensityMeasurements))=[];                    
                    
                    [cellGaussians{iNumberOfGaussians,iIntensity,1} cellGaussians{iNumberOfGaussians,iIntensity,2} cellGaussians{iNumberOfGaussians,iIntensity,3}] = fit_mix_gaussian(matCurrentIntensityMeasurements, iNumberOfGaussians);
                    matIntensityThresholds(iNumberOfGaussians,iIntensity) = cellGaussians{iNumberOfGaussians,iIntensity,1}(1) + (intSDsForThreshold * cellGaussians{iNumberOfGaussians,iIntensity,2}(1));
                end %iIntensity for loop
            end %iNumberOfGaussians for loop
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SAVING AND EXPORTING DATA %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        handles.Measurements.(strObjectName).VirusScreenGaussians{intSetBeingAnalyzed} = cellGaussians;
        handles.Measurements.(strObjectName).VirusScreenThresholds{intSetBeingAnalyzed} = matIntensityThresholds;        

    catch [ErrorMessage, ErrorMessage2] = lasterr;
        error(['Programmo Alert! Check this out: ', ErrorMessage, ErrorMessage2])
    end 

end