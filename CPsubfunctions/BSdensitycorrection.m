function handles = BSdensitycorrection(handles, strObjectName, strCtrlWells)

% load('U:\Data\Users\Berend\density_test\DefaultOUT1.mat');
% strObjectName = 'Nuclei';
% strCtrlWells = 'B05,D06,E07';

%%%%%%%%%%%%%%%%%%%%%
%%% GRID FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%

matRows = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'};
matCols = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};

wellcounter = 0;
str2match = cell(96,1);

xStep = 1440 / 30; %(48)
yStep = 936 / 18; %(52)

% xStep = 720 / 30;
% yStep = 468 / 18;

totalgrid = [0];
infectedgrid = [0];

ctrlwelltotalgrid = [0];
ctrlwellinfectedgrid = [0];
ctrlwellinfectionindexgrid = [0];

%convert ImageNames to something we can index
cellFileNames = {};
for l = 1:length(handles.Measurements.Image.FileNames)
    cellFileNames{l} = char(handles.Measurements.Image.FileNames{l}(1));
end

cellstrCtrlWellsToDo = strread(strCtrlWells,'%s','delimiter',',')';

%first run, just one row, the column loop will loop over the amount of
%control wells, second run, do your average 384 well plate setup.
rowstodo = {[1:1];[1:16]};
%assuming standard 96 well plate 50K setup!
colstodo = {[1:length(cellstrCtrlWellsToDo)];[1:24]};
                    
cellPlateData = {};
cellPlateDataFeatures = {};

intOutOfFocusFeatureNumber = 0;
intOutOfFocusFeatureNumber = strmatch('OutOfFocus', handles.Measurements.(strObjectName).VirusInfectionFeatures);
intDiscardedControlImages = 0;
intTotalControlImages = 0;


for runnumber = 1:2
    wellcounter = 0;
    for rowNum = rowstodo{runnumber} %2:7
        if runnumber == 2
%            disp(['row ', num2str(rowNum),'...'])
        end
        for colNum = colstodo{runnumber} %2:11
            wellcounter = wellcounter + 1;

            if runnumber == 1 %control well analysis
                str2match{wellcounter,1} = strcat(cellstrCtrlWellsToDo{wellcounter}, 'f');
                str2match{wellcounter,2} = strcat(matRows(rowNum), matCols(colNum));                
%                disp(['ctrl well ', cellstrCtrlWellsToDo{wellcounter},'...'])
            else %for loop over all wells (96 well format)
                str2match{wellcounter,1} = strcat(matRows(rowNum), matCols(colNum), 'f');
                str2match{wellcounter,2} = strcat(matRows(rowNum), matCols(colNum));
%                disp(['     col ', num2str(colNum),'...'])                            
            end
            
            FileNameMatches = strfind(cellFileNames, char(str2match{wellcounter,1}));

            if length(find(~cellfun('isempty',FileNameMatches))) > 0
                for k = find(~cellfun('isempty',FileNameMatches))
%                    disp(['         img ', num2str(k),'...'])
                    InfectedNuclei = handles.Measurements.(strObjectName).VirusScreenInfectedNuclei{k};
                    NucleiLocations = handles.Measurements.(strObjectName).Location{k};

                    if length(NucleiLocations) > 0
                        counter = 0;
                        gridindexes = {};
                        gridtotals1d = {};
                        gridbinindexes = {};
                        gridbinnucleinumbers = {};
                        totalinfectedgridbin = 0;

                        for xPos = 0:xStep:1440-1 %-1 to skip last class
                            countery = 0;
                            for yPos = 0:yStep:936-1 %-1 to skip last class
                                counter = counter + 1;
                                gridindexes{counter,1} = find((NucleiLocations(:,1) > xPos) & (NucleiLocations(:,1) < xPos+xStep) & (NucleiLocations(:,2) > yPos) & (NucleiLocations(:,2) < yPos+yStep));
                                gridtotals1d{counter,1} = length(find((NucleiLocations(:,1) > xPos) & (NucleiLocations(:,1) < xPos+xStep) & (NucleiLocations(:,2) > yPos) & (NucleiLocations(:,2) < yPos+yStep)));                
                            end
                        end

                        for bin = 0:30
                            gridbinindexes{bin+1,1} = find(cell2mat(gridtotals1d) == bin);
                            gridbinnucleinumbers{bin+1,1} = cell2mat(gridindexes(gridbinindexes{bin+1,1},1));

                            for i = 1:length(gridbinnucleinumbers{bin+1,1})
                                if find(InfectedNuclei == gridbinnucleinumbers{bin+1,1}(i,1))
                                    totalinfectedgridbin = totalinfectedgridbin + 1;
                                end
                            end

                            totalgrid(bin+1,k) = length(gridbinnucleinumbers{bin+1,1});
                            infectedgrid(bin+1,k) = totalinfectedgridbin;
                            totalinfectedgridbin = 0;
                        end
                    end %for bin
                    
                    
                    if runnumber == 1
                        %determine total number and number of out of focus control images
                        intTotalControlImages = intTotalControlImages + 1;                    
                        if handles.Measurements.(strObjectName).VirusInfection{k}(1,intOutOfFocusFeatureNumber) > 1.3
                            intDiscardedControlImages = intDiscardedControlImages + 1;
                        end
                    end
                    
                end %for k

                
                welltotalgrid{wellcounter} = sum(totalgrid,2);            
                wellinfectedgrid{wellcounter} = sum(infectedgrid,2);
                wellinfectionindexgrid{wellcounter} = sum(infectedgrid,2) ./ sum(totalgrid,2);

                if runnumber == 1
                    %calculate total ctrl correction indexes per bin
                    ctrlwelltotalgrid = ctrlwelltotalgrid + welltotalgrid{wellcounter};         %LD total distribution
                    ctrlwellinfectedgrid = ctrlwellinfectedgrid + wellinfectedgrid{wellcounter};%LD infected distribution
                    ctrlwellinfectionindexgrid = ctrlwellinfectedgrid ./ ctrlwelltotalgrid;     %LD infectionindex distribution

                    intCtrlII = sum(ctrlwellinfectedgrid) / sum(ctrlwelltotalgrid);     % ctrl average infection index
                    intCtrlCN = sum(ctrlwelltotalgrid) / length(cellstrCtrlWellsToDo);  % ctrl average total
                else
                    %apply per bin correction index to calculate expected
                    %infected cells per well
                    temp = (ctrlwellinfectionindexgrid .* welltotalgrid{wellcounter});
                    wellexpectedinfectedgrid(wellcounter) = sum(temp(find(~isnan(temp))));
                    
                    cellPlateData{wellcounter,1} = char(cellFileNames(k));    % file name
                    cellPlateData{wellcounter,2} = char(matRows(rowNum));     % row
                    cellPlateData{wellcounter,3} = char(matCols(colNum));     % column
                    
                    cellPlateData{wellcounter,4} = sum(welltotalgrid{wellcounter});     % total
                    cellPlateData{wellcounter,5} = sum(wellinfectedgrid{wellcounter});  % infected
                    cellPlateData{wellcounter,6} = sum(wellinfectedgrid{wellcounter})/sum(welltotalgrid{wellcounter});    % infection index
                    
                    if char(str2match{wellcounter,1}) == char(strcat(cellstrCtrlWellsToDo{1}, 'f'))
                        cellPlateData{wellcounter,7} = intCtrlII;
                        cellPlateData{wellcounter,8} = intCtrlCN;
                    else
                        cellPlateData{wellcounter,7} = [];
                        cellPlateData{wellcounter,8} = [];                       
                    end

                    cellPlateData{wellcounter,9} = cellPlateData{wellcounter,6}/intCtrlII; % relative infection index
                    cellPlateData{wellcounter,10} = cellPlateData{wellcounter,4}/intCtrlCN;% relative cell number

                    cellPlateData{wellcounter,11} = log2(cellPlateData{wellcounter,9});  % log2(relative infection index)
                    cellPlateData{wellcounter,12} = log2(cellPlateData{wellcounter,10}); % log2(relative cell number)

                    cellPlateData{wellcounter,13} = wellexpectedinfectedgrid(wellcounter);                              % ctrl LDcorr expected infected
                    cellPlateData{wellcounter,14} = wellexpectedinfectedgrid(wellcounter)/cellPlateData{wellcounter,4};   % ctrl LDcorr infection index
                    cellPlateData{wellcounter,15} = cellPlateData{wellcounter,6}/cellPlateData{wellcounter,14};   % LDcorr relative infection index                    
                    cellPlateData{wellcounter,16} = log2(cellPlateData{wellcounter,15});   % ctrl LDcorr infection index                                        
                end                
            end %if length(matches) > 0
            infectedgrid = [0];
            totalgrid = [0];
            temp = [0];
        end
    end
end

if not(isempty(cellPlateData))
    matFullIndices = find(~cellfun('isempty',cellPlateData(:,1)));
    handles.Measurements.(strObjectName).VirusScreenLDCorrectedPlateData = cellPlateData(matFullIndices,:);
else
    handles.Measurements.(strObjectName).VirusScreenLDCorrectedPlateData = [];
end
cellPlateDataFeatures = {'File','Row','Column','Total','Infected','II','CtrlII','CtrlCN','RII','RCN','log2(RII)','log2(RCN)','ctrlLDCorrInfected','ctrlLDCorrII','RLDCorrII','log2(RLDCorrII)'};
handles.Measurements.(strObjectName).VirusScreenLDCorrectedPlateDataFeatures = cellPlateDataFeatures;
handles.Measurements.(strObjectName).VirusScreenDiscardedControls{1}(1,1) = intDiscardedControlImages ;
handles.Measurements.(strObjectName).VirusScreenTotalControls{1}(1,1) = intTotalControlImages;
