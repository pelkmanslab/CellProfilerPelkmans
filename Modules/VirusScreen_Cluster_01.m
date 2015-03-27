function handles = VirusScreen_Cluster_01(handles)

warning off all

% Help for the VirusScreen_Cluster_01 module:
% Category: Other
%
% SHORT DESCRIPTION:
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

%%%VariableRevisionNumber = 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PREPARATION OF VARIABLES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    intSetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

    if (intSetBeingAnalyzed == 1)
        handles.Settings.VariableNames{CurrentModuleNum} = {
            'strObjectName',...
            'strImageName',...        
        };
    end

    strIntensityMeasurementList = '1,2,3,4,5,6,7,8,9,10,11';
    matIntensityMeasurements = str2num(strrep(strIntensityMeasurementList, ',', ' '));
    
    intSDsForThreshold = 4;

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
                    [cellGaussians{iNumberOfGaussians,iIntensity,1} cellGaussians{iNumberOfGaussians,iIntensity,2} cellGaussians{iNumberOfGaussians,iIntensity,3}] = fit_mix_gaussian(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)), iNumberOfGaussians);
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
