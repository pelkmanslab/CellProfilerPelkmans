function handles = VirusScreen_Auto_07_sensitive(handles)

warning off all

% Help for the VirusScreen_Auto_07_sensitive module:
% Category: Other
%
% SHORT DESCRIPTION:
% VirusScreen_Auto_07_sensitive Does double gaussian fitting for two measurements, selects
% those nuclei that are hits for both groups and only scores those. It furthermore
% does some postanalysis, and calls to ExportVirusScreen_01.m and
% ExportVirusScreenLogfile.m for data export.
% SCORING 3 SD'S FROM BACKGROUND
% *************************************************************************
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Which object would you like to use for the infection measurement?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Which image contains the virus signal measurements?
%infotypeVAR02 = imagegroup
%inputtypeVAR02 = popupmenu
strImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Are your no-virus control images infected?
%choiceVAR03 = No
%choiceVAR03 = Yes
%inputtypeVAR03 = popupmenu
strInfectedControls = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Do you want to apply local density correction?
%choiceVAR04 = Yes
%choiceVAR04 = No
strLocalDensityCorrection = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu    

%textVAR05 = Which are your local density control wells (comma separated)?
%defaultVAR05 = B05,B06,B07
strLocalDensityCtrlWells = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%%%VariableRevisionNumber = 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PREPARATION OF VARIABLES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = handles.Current.SetBeingAnalyzed;

    if (k == 1)
        handles.Settings.VariableNames{CurrentModuleNum} = {
            'strObjectName',...
            'strImageName',...
            'strInfectedControls'...
            'strLocalDensityCorrection'...
            'strLocalDensityCtrlWells'...            
        };
    end

    intOutOfFocusCutoff = 1.3;

    strIntensityMeasurementList = '2,6';
    intSDsForThreshold = 3;

    strIntensityFieldName = ['Intensity_', strImageName];

    % Check if the module LoadImagesControlFirst is loaded, and check for the
    % control images
    intLoadImagesControlFirstModuleNumber = strmatch('LoadImagesControlFirst', handles.Settings.ModuleNames);
    strControlImagesRegExp = char(handles.Settings.VariableValues{intLoadImagesControlFirstModuleNumber,15});

    strCurrentImageName = handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{2};

    intTotalCells = handles.Measurements.Image.ObjectCount{k}(1,1);

    
    intUpper90Diameter = [];
    intMaxDiameter = [];
    intObjectCoverage = [];
    intDiameterExcludedObjects = [];
    intBorderObjects = [];
    
    try
        intUpper90Diameter = handles.Measurements.Image.ObjectDiameter{k}(1,1);
        intMaxDiameter = handles.Measurements.Image.ObjectDiameter{k}(1,2);
        intObjectCoverage = handles.Measurements.Image.ObjectDiameter{k}(1,3);
        intDiameterExcludedObjects = handles.Measurements.Image.ObjectDiameter{k}(1,4);
        intBorderObjects = handles.Measurements.Image.ObjectDiameter{k}(1,5);
    end

    IntensityMeasurements = handles.Measurements.(strObjectName).(strIntensityFieldName){k}; %(:,intFirstIntensityMeasurement);

    intFlagOutOfFocus = 0;
    if intUpper90Diameter > intMaxDiameter
        intFlagOutOfFocus = (intUpper90Diameter / intMaxDiameter);
    end

    boolDiscardImage = 0;
    boolControlImage = 0;
    intNumberOfGaussians = 2;
    intTotalInfectedCells1 = 0;
    intInfectionIndex = 0;
    doubleHits1 = [];
    cellSingleHits = {};
    matIntensityMeasurements = str2num(strrep(strIntensityMeasurementList, ',', ' '));
    matIntensityThresholds = [];
    matPreliminaryHits = [];
    strThresholdUsed = cell(1,length(matIntensityMeasurements));
    strFigureTitle = '';
    cellGaussians = cell(length(matIntensityMeasurements), 3);
    cellIncludedThresholds = cell(1,length(matIntensityMeasurements));

    if ~isfield(handles.Measurements.(strObjectName),'VirusScreenIncludedThresholds')
        handles.Measurements.(strObjectName).VirusScreenIncludedThresholds = {};
    else
        cellIncludedThresholds = handles.Measurements.(strObjectName).VirusScreenIncludedThresholds;
    end

    if ~isfield(handles.Measurements.(strObjectName),'VirusScreenThresholds')
        handles.Measurements.(strObjectName).VirusScreenThresholds = {};
    end

    if ~isfield(handles.Measurements.(strObjectName),'VirusScreenWindowCheck')
        handles.Measurements.(strObjectName).VirusScreenWindowCheck = {[0]};
    end

    intWindowOpenForNumberOfCycles = handles.Measurements.(strObjectName).VirusScreenWindowCheck{1}(1,1);

%     try

        if intTotalCells > 1        

            if not(isempty(regexp(strCurrentImageName,strControlImagesRegExp,'ONCE')))
                boolControlImage = 1;
                if strcmp(strInfectedControls,'Yes')
                    intNumberOfGaussians = 2;
                else
                    intNumberOfGaussians = 1;                
                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% INFECTIVITY CALCULATION %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % Do all two gaussian fittings
            for iIntensity = 1:length(matIntensityMeasurements)
                [cellGaussians{iIntensity,1} cellGaussians{iIntensity,2} cellGaussians{iIntensity,3}] = fit_mix_gaussian(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)), intNumberOfGaussians);
                matIntensityThresholds(1,iIntensity) = cellGaussians{iIntensity,1}(1) + (intSDsForThreshold * cellGaussians{iIntensity,2}(1));
                if isnan(matIntensityThresholds(1,iIntensity))
                    if length(cellIncludedThresholds{iIntensity}) > 0
                        matPreliminaryHits{iIntensity} = find(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)) > mean(cellIncludedThresholds{iIntensity}));
                        strThresholdUsed{1,iIntensity} = 'average, because current failed';
                    else
                        matPreliminaryHits{iIntensity} = [];
                        strThresholdUsed{1,iIntensity} = 'none, because current failed and no average yet';
                    end
                else
                    if length(cellIncludedThresholds{iIntensity}) > 0
                        %if there is an average Threshold
                        if (matIntensityThresholds(1,iIntensity) > (1.5 * mean(cellIncludedThresholds{iIntensity})))       
                            %and if the current Threshold is 1,5x bigger than
                            %the mean, use the mean...
                            matPreliminaryHits{iIntensity} = find(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)) > mean(cellIncludedThresholds{iIntensity}));
                            strThresholdUsed{1,iIntensity} = 'average, because current is too high';
                        elseif (matIntensityThresholds(1,iIntensity) < min(cellIncludedThresholds{iIntensity}))
                            %and if the current Threshold is smaller then
                            %the lowest control Threshold, use the mean...
                            matPreliminaryHits{iIntensity} = find(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)) > mean(cellIncludedThresholds{iIntensity}));
                            strThresholdUsed{1,iIntensity} = 'average, because current is too low';                        
                        else
                            % otherwise, use current Threshold.
                            matPreliminaryHits{iIntensity} = find(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)) > matIntensityThresholds(1,iIntensity));
                            strThresholdUsed{1,iIntensity} = 'current fit';
                        end
                    else
                        %always use current Threshold if there is no average
                        %Threshold yet.
                        matPreliminaryHits{iIntensity} = find(IntensityMeasurements(:,matIntensityMeasurements(1,iIntensity)) > matIntensityThresholds(1,iIntensity));
                        strThresholdUsed{1,iIntensity} = 'current fit';
                    end
                end

                if boolControlImage && (intTotalCells > 30) && not(isnan(matIntensityThresholds(1,iIntensity)))
                    matIntensityThresholds(2,iIntensity) = 1;
                    cellIncludedThresholds{iIntensity}(end+1) = matIntensityThresholds(1,iIntensity);
                else
                    matIntensityThresholds(2,iIntensity) = 0;
                end           

            end %iIntensity for loop

            for ix = 1:intTotalCells
                arrScoreList(ix,1) = 0;
                arrScoreList(ix,2) = 0;    

                for iy = 1:length(matPreliminaryHits{1})
                    if matPreliminaryHits{1}(iy,1) == ix
                        arrScoreList(ix,1) = 1;
                    end
                end
                for iz = 1:length(matPreliminaryHits{2})     
                    if matPreliminaryHits{2}(iz,1) == ix
                        arrScoreList(ix,2) = 1;
                    end
                end
                arrScoreList(ix,5) = ix;
            end

            arrScoreList(:,3) = handles.Measurements.(strObjectName).Location{k}(:,1);
            arrScoreList(:,4) = handles.Measurements.(strObjectName).Location{k}(:,2);

            for ib = 1:length(matIntensityMeasurements)
                cellSingleHits{ib} = arrScoreList(find(arrScoreList(:,ib) == 1), :);
            end

            doubleHits1 = arrScoreList(find(arrScoreList(:,1) == 1 & arrScoreList(:,2) == 1), :);
            intTotalInfectedCells1 = length(doubleHits1(:,1));

            intInfectionIndex = (intTotalInfectedCells1 / intTotalCells);

            if (intFlagOutOfFocus > intOutOfFocusCutoff) %|| ((intTotalCells < intMinNucleusCutoff) && (sboolPossibleSinglePopulation == 1 || fboolPossibleSinglePopulation == 1))
                boolDiscardImage = 1;
            end

        end


        %%%%%%%%%%%%%%%%%%%%%%%
        %%% DISPLAY RESULTS %%%
        %%%%%%%%%%%%%%%%%%%%%%%

        ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

        if intWindowOpenForNumberOfCycles == 100
            if any(findobj == ThisModuleFigureNumber)
                close(ThisModuleFigureNumber);
            end
        end    

        if any(findobj == ThisModuleFigureNumber)
            intWindowOpenForNumberOfCycles = intWindowOpenForNumberOfCycles + 1;
            CPfigure(handles,'Image',ThisModuleFigureNumber);
            if intTotalCells > 1

                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Display all plots %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%        

                % Figure Title
                if boolControlImage
                    strFigureTitle = [strCurrentImageName, '  -  ', num2str(k),'/',num2str(handles.Current.NumberOfImageSets), '  (CONTROL)'];
                else
                    strFigureTitle = [strCurrentImageName, '  -  ', num2str(k),'/',num2str(handles.Current.NumberOfImageSets)];
                end
                if intWindowOpenForNumberOfCycles > 2
                    uicontrol('style','text','fontsize',handles.Preferences.FontSize,'fontweight','bold','HorizontalAlignment','center','String',['ATTENTION: THIS WINDOW WILL CLOSE AUTOMATICALLY SOON'],'position',[10 10 400 20],'BackgroundColor',[.7 .7 .9]);                
                end
                uicontrol('style','text','units','normalized','fontsize',handles.Preferences.FontSize,'fontweight','bold','HorizontalAlignment','center','String',strFigureTitle,'position',[.1 0.95 .8 .05],'BackgroundColor',[.7 .7 .9]);

                subplot(4,4,1);

                for iPlots = 1:length(matIntensityMeasurements)
                    %plot-title
                    strPlotName = [handles.Measurements.(strObjectName).([strIntensityFieldName,'Features']){matIntensityMeasurements(1,iPlots)}]; 
                    if intNumberOfGaussians == 1
                        plot_mix_gaussian2(cellGaussians{iPlots,1},cellGaussians{iPlots,2},[1],IntensityMeasurements(:,matIntensityMeasurements(1,iPlots)),4,4,iPlots,strPlotName,[matIntensityThresholds(1,iPlots), mean(cellIncludedThresholds{iPlots}), min(cellIncludedThresholds{iPlots})]);
                    else
                        plot_mix_gaussian2(cellGaussians{iPlots,1},cellGaussians{iPlots,2},[0.5, 0.5],IntensityMeasurements(:,matIntensityMeasurements(1,iPlots)),4,4,iPlots,strPlotName,[matIntensityThresholds(1,iPlots), mean(cellIncludedThresholds{iPlots}), min(cellIncludedThresholds{iPlots})]);
                    end
                    %plot-subscript
                    uicontrol('style','text','units','normalized','HorizontalAlignment','center','fontsize',7,'String',['(',strThresholdUsed{1,iPlots},') (ctrl_n = ', num2str(length(cellIncludedThresholds{iPlots})),')'],'position',[-.1+(iPlots*.21) 0.72 .2 .015],'BackgroundColor',[.7 .7 .9]);
                end

                if isfield(handles.Pipeline,'SortedMaxValuesMatrixOrigGreen')
                    try
                        %plot intensity histogram if this data is present... :)
                        %strImageName... hardcoded OrigGreen... man I'm just getting
                        %lasier and lasier... :)
                        fieldname = ['SortedMaxValuesMatrix', 'OrigGreen'];
                        matIntensityValues = handles.Pipeline.(fieldname);
                        fieldname = ['MaxPixelValue', 'OrigGreen'];
                        HighestPixelOrig = handles.Pipeline.(fieldname);
                        iPlots = length(matIntensityMeasurements)+1;            
                        strPlotName = '(Max Intensities of All Images)'; 
                        subplot(4,4,iPlots);
                        plot(matIntensityValues);
                        title(strPlotName, 'fontsize', 7,'Color',[0.3, 0.3, 0.3]);
                        vline((0.3*length(matIntensityValues)),'k','');
                        hline(HighestPixelOrig,'k','');            
                        uicontrol('style','text','units','normalized','HorizontalAlignment','center','fontsize',7,'String',['(rescaled to ',num2str(HighestPixelOrig),')'],'position',[-.1+(iPlots*.21) 0.72 .2 .015],'BackgroundColor',[.7 .7 .9]);
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Display Virus image %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(4,4,[5:16]);                     
                OrigImage = CPretrieveimage(handles,strImageName,ModuleName,'MustBeGray','CheckScale');
                %[BS] Removed these guys to try to improve display
                %issues... Like overwriting other figure windows...
                %
                %ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
                %CPfigure(handles,'Image',ThisModuleFigureNumber);
                CPimagesc(OrigImage,handles);
                colormap(gray);
                text(doubleHits1(:,3) , doubleHits1(:,4), cellstr(num2str(doubleHits1(:,5))),...
                    'HorizontalAlignment','center', 'color', [0 .8 0],'fontsize',7);         

                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Display text stuff %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%

                TextString = {['TotalCells = ',num2str(double(intTotalCells))], ['InfectedCells = ', num2str(double(intTotalInfectedCells1))], ['InfectionIndex = ', num2str(double(intInfectionIndex))],[' '],  ['FlagOutOfFocus = ', num2str(double(intFlagOutOfFocus))], ['Discarded ', strObjectName, ' Objects = ',num2str(double(intDiameterExcludedObjects + intBorderObjects))], ['ObjectCoverage = ',num2str(double(intObjectCoverage))]};
                for i = 1:length(TextString)
                    uicontrol('style','text','units','normalized','fontsize',8,'HorizontalAlignment','left','String',TextString{i},'position',[.75 0.92-0.02*i .2 .02],'BackgroundColor',[.7 .7 .9]);
                end

            end
        else %any(findobj == ThisModuleFigureNumber)
            intWindowOpenForNumberOfCycles = 0;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SAVING AND EXPORTING DATA %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        drawnow

        handles.Measurements.(strObjectName).VirusScreenWindowCheck{1}(1,1) = intWindowOpenForNumberOfCycles;    
        handles.Measurements.(strObjectName).VirusScreenThresholds{k} = matIntensityThresholds;
        handles.Measurements.(strObjectName).VirusScreenIncludedThresholds = cellIncludedThresholds;
        if not(isempty(doubleHits1))
            handles.Measurements.(strObjectName).VirusScreenInfectedNuclei{k} = doubleHits1(:,5);
        else
            handles.Measurements.(strObjectName).VirusScreenInfectedNuclei{k} = [];            
        end

        handles = CPaddmeasurements(handles,strObjectName,'VirusInfection','TotalCells',intTotalCells);
        handles = CPaddmeasurements(handles,strObjectName,'VirusInfection','TotalInfectedCells',intTotalInfectedCells1);
        handles = CPaddmeasurements(handles,strObjectName,'VirusInfection','OutOfFocus',intFlagOutOfFocus);
        handles = CPaddmeasurements(handles,strObjectName,'VirusInfection','DiscardImage',boolDiscardImage);        

        if (k == handles.Current.NumberOfImageSets) && strcmp(strLocalDensityCorrection, 'Yes')
            handles = BSdensitycorrection(handles, strObjectName, strLocalDensityCtrlWells);
        end
% 
%     catch [ErrorMessage, ErrorMessage2] = lasterr;
%         error(['Programmo Alert! Check this out: ', ErrorMessage, ErrorMessage2])
%     end 

    if (k == handles.Current.NumberOfImageSets)
        % on last cycle, produce export files.
%        if  strcmp(strLocalDensityCorrection, 'Yes')
%            ExportVirusScreenLocalDensity_01(handles);
%        end
        ExportVirusScreen_01(handles);
        ExportVirusScreenLogfile_02(handles);
    end
end
