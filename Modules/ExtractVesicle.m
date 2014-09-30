function handles = ExtractVesicle(handles)

% Help for the ExtractVesicle module:
% Category: Other
%
% SHORT DESCRIPTION:
% Identifies objects (e.g. cell edges) using "seed" objects identified by
% an Identify Primary module (e.g. nuclei).
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4904 $

% 22.07.2008 Pekka Ruusuvuori
% 07.11.2008 PR - results for multiple channels enabled
% 16.02.2009 Sharif Chowdhury - Riples K function added and output is saved
% 17.02.2009 Sharif Chowdhury - Some Presentation of measurement is changed
% 18.02.2009 Sharif Chowdhury - Output is produced according to the object
% names

% using CPaddmeasurements function & detection module is removed from here
% compatible CP version: 1.0.4553

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object Name?
%defaultVAR01 = Vesicles
%infotypeVAR01 = objectgroup indep

%textVAR02 = Segmented Secondary Outline Name?
%defaultVAR02 = Cells
%infotypeVAR02 = objectgroup indep

%textVAR03 = Which image contains the signal?
%infotypeVAR03 = imagegroup
%inputtypeVAR03 = popupmenu

%textVAR04 = Riple's Starting Radius ?
%defaultVAR04 = 1
%infotypeVAR04 = objectgroup indep

%textVAR05 = Riple's Ending Radius (zero for disabling the feature) ?
%defaultVAR05 = 0
%infotypeVAR05 = objectgroup indep

%%%VariableRevisionNumber = 1




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
drawnow

%MainModNo = 1;


strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
strMaskObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
strImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
rad1 = str2num( char( handles.Settings.VariableValues{CurrentModuleNum,4}) );
if rad1<1
    rad1 = 1;
end
rad2 = str2num( char( handles.Settings.VariableValues{CurrentModuleNum,5}) );



%R = handles.Pipeline.(char(handles.Settings.VariableValues{MainModNo,7}));
%G = handles.Pipeline.(char(handles.Settings.VariableValues{MainModNo,5}));
%B = handles.Pipeline.(char(handles.Settings.VariableValues{MainModNo,3}));

%R = handles.Pipeline.OrigRed;
%G = handles.Pipeline.OrigGreen;
%B = handles.Pipeline.OrigBlue;
%I = uint16(cat(3,R,G,B));
%I = imread(strcat(path,imnames(imind,:),num2str(channel),'.tif'));
%L = imread([path imnames(imind,:) '0_SegmentedCells.tif']);
fieldname = ['Segmented',strMaskObjectName];
L = uint16(CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale'));
G = handles.Pipeline.(strImageName);

%figure(imind)
%h(1) = subplot(221)
%imshow(L,[]);
%h(2) = subplot(222)
%imshow(log(double(I))/max(max(max(log(double(I))))),[]);


imind = handles.Current.SetBeingAnalyzed;
% channel = 1;
% stats = regionprops(L, 'PixelIdxList', 'PixelList');
% b = bwboundaries(L);
% centerpoints = handles.Measurements.Cells.Location{imind};

% vesicle detection
th_sensitivity = 1.6; % around 1.5 seems to be OK 
visualize = 0;
%%[ves,CF,vesres] = brightspots(G,th_sensitivity,visualize);

fieldname = ['Segmented', strObjectName];
ves= CPretrieveimage(handles,fieldname,ModuleName,'DontCheckColor','DontCheckScale',size(L));

%figure(imind)
%h(3) = subplot(223)
%imshow(ves)
%linkaxes(h) 

% extract cell level features
[features,allvesicledata,fnames, vfnames, riplyMatrix] = extract_vesicle_features( double(L), G, ves, rad1, rad2);


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

% if any(findobj == ThisModuleFigureNumber)
%     %%% Activates the appropriate figure window.
%     CPfigure(handles,'Image',ThisModuleFigureNumber);
%      CPimagesc(vesres,handles);
%     title(['Detected vesicles, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

featurefieldname = [ModuleName '_',strImageName,'Features'];
fieldname = [ModuleName '_',strImageName];
%handles.Measurements.Image.(fieldname)(handles.Current.SetBeingAnalyzed) = {[TotalIntensity MeanIntensity TotalArea]};

%handles.Measurements.(objectname).(modulename_RescaledRed){imind}(n_objects, m_measurements)
%handles.Measurements.(objectname).(modulename_RescaledRedFeatures) = {'column1_description','column2_description'};
%handles.Measurements.(strObjectName).(fieldname){imind} = features;
%handles.Measurements.(strObjectName).(featurefieldname) = fnames;


% % % featurefieldname = ['Intensity_',ImageName,'Features'];
% % % 
% % % fieldname = ['Intensity_',ImageName];
% % % 
% % % handles.Measurements.Image.(featurefieldname) = 
% % % {'TotalIntensity','MeanIntensity','TotalArea'};
% % % 
% % % handles.Measurements.Image.(fieldname)(handles.Current.SetBeingAnalyzed) = 
% % % {[TotalIntensity MeanIntensity TotalArea]};


%%function handles = CPaddmeasurements(handles,Object,Measure,Feature,Data)
%%FeaturesField = [Measure,'Features'];
% % handles.Measurements.(Object).(FeaturesField) = {Feature};
% % handles.Measurements.(Object).(Measure){handles.Current.SetBeingAnalyzed} = Data;

 if  isfield(handles.Measurements,strMaskObjectName) < 1
     handles.Measurements.(strMaskObjectName) = cell(0);
 end

 handles =CPaddmeasurements(handles,strMaskObjectName,['Number_of_',strObjectName,'_', strImageName],'Total_Count',features(:,1));
 handles =CPaddmeasurements(handles,strMaskObjectName,[strObjectName,'_Size_', strImageName],'Average',features(:,2));
 handles =CPaddmeasurements(handles,strMaskObjectName,[strObjectName,'_Intensity_', strImageName],'Average',features(:,3));
 handles =CPaddmeasurements(handles,strMaskObjectName,['Number_of_',strObjectName,'_Per_Unit_Cell_Area_',strImageName],'Average',features(:,4));
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_of_',strObjectName,'_to_Cell_Outline_', strImageName],'Average',features(:,5));
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_of_',strObjectName,'_to_Cell_Outline_', strImageName], 'STD',features(:,6));
 
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_of_',strObjectName,'_to_Cell_Centre_',strImageName],'Average',features(:,7));
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_of_',strObjectName,'_to_Cell_Centre_',strImageName], 'STD',features(:,8));
 
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_Between_',strObjectName,'_', strImageName],'Average',features(:,9));
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_Between_',strObjectName,'_', strImageName],'SRD',features(:,10));
 handles =CPaddmeasurements(handles,strMaskObjectName,['Distance_Between_',strObjectName,'_', strImageName], 'Minimum',features(:,11));

 if rad2 > 0 && rad2>=rad1
     for i=rad1:rad2 
        handles =CPaddmeasurements(handles,strMaskObjectName,['RiplesK_',strObjectName,'_', strImageName], ['Radius -' num2str(i) ],riplyMatrix(:,i));
     end
 end
 
 
 
 if  isfield(handles.Measurements,strObjectName) < 1
     handles.Measurements.(strObjectName) = cell(0);
 end
     
 handles =CPaddmeasurements(handles, strObjectName,['Distance_', strImageName], 'To_Cell_Outline',allvesicledata(:,1));
 handles =CPaddmeasurements(handles, strObjectName,['Distance_', strImageName], 'To_Cell_Centre',allvesicledata(:,2));
 handles =CPaddmeasurements(handles, strObjectName,['Intensity_', strImageName], 'Mean',allvesicledata(:,3));
 handles =CPaddmeasurements(handles, strObjectName,['Area_', strImageName], 'Total',allvesicledata(:,4));
 handles =CPaddmeasurements(handles, strObjectName,['Orientation_', strImageName], 'Angle',allvesicledata(:,5));
 handles =CPaddmeasurements(handles, strObjectName,['Parent_Cell_', strImageName], 'CellID',allvesicledata(:,6));
 

end