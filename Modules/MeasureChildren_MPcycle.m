function handles = MeasureChildren_MPcycle(handles)

% Help for the Measure Texture module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Pools measurements of children to create measurments for each parent. The
% measurments include mean/median and var as well as higher central
% moments. Values with NaNs are excluded from the analysis.
% *************************************************************************
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles); %#ok<NASGU>

%textVAR01 = What objects are the children objects (subobjects)?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
SubObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What are the parent objects?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
ParentName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What feature do you want to use (option: 'all')?
%defaultVAR03 = all
FeatureName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Do you want to rename the children objects?
%choiceVAR04 = Do not use
%infotypeVAR04 = objectgroup
%inputtypeVAR04 = popupmenu
NewObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = What did you call the images the objects are related to?
%infotypeVAR05 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu


%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drawnow

if strcmp(FeatureName,'all')
    AllFeatures = fieldnames(handles.Measurements.(SubObjectName));
    AllFeatures = AllFeatures(cellfun(@(x) isempty(strfind(x,'Features')),AllFeatures)); % this could be done better
    AllFeatures = AllFeatures(cellfun('isempty',regexp(AllFeatures,'Parent|Location|Localization')));
else
    AllFeatures = FeatureName;
end

for i = 1:length(AllFeatures)
    
    fprintf('%s: pools measurements of children ''%s'' for parent ''%s'': feature ''%s''\n',mfilename,SubObjectName,ParentName,AllFeatures{i})
    
    % Import Features of Interst of Child
    matImpSpotFeatures =  handles.Measurements.(SubObjectName).(AllFeatures{i}){handles.Current.SetBeingAnalyzed};
    matImpSpotFeatureDescription =  handles.Measurements.(SubObjectName).([AllFeatures{i} 'Features']);
    
    % Get number of parent objects
    column = find(strcmpi(handles.Measurements.Image.ObjectCountFeatures,ParentName),1,'first');
    matTempParentCount =  handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed};
    matParentCount = matTempParentCount(:,column); clear matTempParentLocation;
    
    % Get the ID of the parent for each child
    matParentObjectAll = handles.Measurements.(SubObjectName).Parent{handles.Current.SetBeingAnalyzed};
    columnParentOfInterst = find(strcmpi(handles.Measurements.(SubObjectName).ParentFeatures,ParentName),1,'first');
    matParentObject = matParentObjectAll(:,columnParentOfInterst);
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %%% DATA ANALYSIS %%%
    %%%%%%%%%%%%%%%%%%%%%
    drawnow
    
    numParents = max([max(matParentObject(:)), matParentCount]);
    
    % initialize output
    PCNaNMean = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNMedian = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNVar = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNMoment3 = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNMoment4 = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNMoment5 = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNMoment6 = nan(numParents,size(matImpSpotFeatures,2));
    PCNaNStd = nan(numParents,size(matImpSpotFeatures,2));
    
    
    if any(matParentObject>0) && (numParents>0);
        
        % sort measurments of children according to parents. this will speed up
        % the later calculation
        [smatParentObject sIX] = sort(matParentObject); % sort children according to their parents
        
        smatImpSpotFeatures = matImpSpotFeatures(sIX,:); % sort data of children so that data of siblings is next to each other. This allows to use a fast strategy for filtering
        
        [usmatParentObject bFirst] = unique(smatParentObject,'first'); % get ID of first child
        [~, bLast] = unique(smatParentObject,'last'); % last child
        
        for j=1:numParents % for each parent
            f = usmatParentObject == j; % identify children
            if any(f)
                CurrData = smatImpSpotFeatures(bFirst(f):bLast(f),:); % import data of all children
                
                PCNaNMean(j,:) =        nanmean(CurrData,1);
                PCNaNMedian(j,:) =      nanmedian(CurrData,1);
                PCNaNVar(j,:) =         nanvar(CurrData,[],1);
                PCNaNStd(j,:) =         nanstd(CurrData,[],1);
                
                
                for k=1:size(smatImpSpotFeatures,2) % loop through each feature. note that it is possible that only some features have nans.
                    
                    ff = ~(isnan(CurrData(:,k)));
                    CurrColumnData = CurrData(ff,k);
                    
                    PCNaNMoment3(j,k) =     moment(CurrColumnData,3,1);
                    PCNaNMoment4(j,k) =     moment(CurrColumnData,4,1);
                    PCNaNMoment5(j,k) =     moment(CurrColumnData,5,1);
                    PCNaNMoment6(j,k) =     moment(CurrColumnData,6,1);
                    
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%
    %%% SAVE RESULTS %%%
    %%%%%%%%%%%%%%%%%%%%
    
    % remove image name tag (we add it later and want to prevent it to
    % occur twice in the string)
    AllFeatures{i} = strrep(AllFeatures{i}, sprintf('_%s', ImageName), '');
    
    % note that for each parent derived measuremnt is saved in new file.
    
    nameMeasurementNaNMean=     [AllFeatures{i} '_ObjectMean'];
    nameMeasurementNaNMedian =  [AllFeatures{i} '_ObjectMedian'];
    nameMeasurementNaNVar =     [AllFeatures{i} '_ObjectVariance'];
    nameMeasurementNaNStd =     [AllFeatures{i} '_ObjectStd'];
    nameMeasurementNaNMom3 =    [AllFeatures{i} '_Object3rdMoment'];
    nameMeasurementNaNMom4 =    [AllFeatures{i} '_Object4thMoment'];
    nameMeasurementNaNMom5 =    [AllFeatures{i} '_Object5thMoment'];
    nameMeasurementNaNMom6 =    [AllFeatures{i} '_Object6thMoment'];

    

    if strcmp(NewObjectName,'Do not use')
        error('define name of children objects')
    else
        ObjectName = NewObjectName;
    end
    
    if handles.Current.SetBeingAnalyzed==1
        
        matImpSpotFeatureDescriptionNaNMean = cellfun(@(x)  [x 'nameMeasurementNaNMean'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNMedian = cellfun(@(x)[x 'nameMeasurementNaNMedian'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNVar = cellfun(@(x)   [x 'nameMeasurementNaNVar'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNMom3 = cellfun(@(x)  [x 'nameMeasurementNaNMom3'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNMom4 = cellfun(@(x)  [x 'nameMeasurementNaNMom4'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNMom5 = cellfun(@(x)  [x 'nameMeasurementNaNMom5'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNMom6 = cellfun(@(x)  [x 'nameMeasurementNaNMom6'], matImpSpotFeatureDescription,'UniformOutput',false);
        matImpSpotFeatureDescriptionNaNStd = cellfun(@(x)   [x 'nameMeasurementNaNStd'], matImpSpotFeatureDescription,'UniformOutput',false);
        
        
        handles.Measurements.(ObjectName).([nameMeasurementNaNMean, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNMedian, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNVar, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom3, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom4, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom5, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom6, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        handles.Measurements.(ObjectName).([nameMeasurementNaNStd, '_', ImageName]) = cell(1,handles.Current.NumberOfImageSets);
        
        
        handles.Measurements.(ObjectName).([nameMeasurementNaNMean, '_', ImageName,      'Features']) =    matImpSpotFeatureDescriptionNaNMean;
        handles.Measurements.(ObjectName).([nameMeasurementNaNMedian, '_', ImageName,    'Features']) =    matImpSpotFeatureDescriptionNaNMedian;
        handles.Measurements.(ObjectName).([nameMeasurementNaNVar, '_', ImageName,       'Features']) =    matImpSpotFeatureDescriptionNaNVar;
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom3, '_', ImageName,      'Features']) =    matImpSpotFeatureDescriptionNaNMom3;
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom4, '_', ImageName,      'Features']) =    matImpSpotFeatureDescriptionNaNMom4;
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom5, '_', ImageName,      'Features']) =    matImpSpotFeatureDescriptionNaNMom5;
        handles.Measurements.(ObjectName).([nameMeasurementNaNMom6, '_', ImageName,      'Features']) =    matImpSpotFeatureDescriptionNaNMom6;
        handles.Measurements.(ObjectName).([nameMeasurementNaNStd, '_', ImageName,       'Features']) =    matImpSpotFeatureDescriptionNaNStd;
        
        
    end
    
    % Save Measurements
    handles.Measurements.(ObjectName).([nameMeasurementNaNMean, '_', ImageName]){handles.Current.SetBeingAnalyzed} =      PCNaNMean;
    handles.Measurements.(ObjectName).([nameMeasurementNaNMedian, '_', ImageName]){handles.Current.SetBeingAnalyzed} =    PCNaNMedian;
    handles.Measurements.(ObjectName).([nameMeasurementNaNVar, '_', ImageName]){handles.Current.SetBeingAnalyzed} =       PCNaNVar;
    handles.Measurements.(ObjectName).([nameMeasurementNaNMom3, '_', ImageName]){handles.Current.SetBeingAnalyzed} =      PCNaNMoment3;
    handles.Measurements.(ObjectName).([nameMeasurementNaNMom4, '_', ImageName]){handles.Current.SetBeingAnalyzed} =      PCNaNMoment4;
    handles.Measurements.(ObjectName).([nameMeasurementNaNMom5, '_', ImageName]){handles.Current.SetBeingAnalyzed} =      PCNaNMoment5;
    handles.Measurements.(ObjectName).([nameMeasurementNaNMom6, '_', ImageName]){handles.Current.SetBeingAnalyzed} =      PCNaNMoment6;
    handles.Measurements.(ObjectName).([nameMeasurementNaNStd, '_', ImageName]){handles.Current.SetBeingAnalyzed} =       PCNaNStd;
    
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
end
