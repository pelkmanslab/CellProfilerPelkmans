function handles = ExtractCell(handles)

% Help for the ExtractCell module:
% Category: Other
%
% SHORT DESCRIPTION:
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4904 $

% 22.07.2008 Pekka Ruusuvuori
% 07.11.2008 PR - results for multiple channels enabled

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%defaultVAR01 = Cells
%infotypeVAR01 = objectgroup indep
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Which image contains the signal?
%infotypeVAR02 = imagegroup
%inputtypeVAR02 = popupmenu
strImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
drawnow

%MainModNo = 1;

% R = handles.Pipeline.(char(handles.Settings.VariableValues{MainModNo,7}));
% G = handles.Pipeline.(char(handles.Settings.VariableValues{MainModNo,5}));
% B = handles.Pipeline.(char(handles.Settings.VariableValues{MainModNo,3}));
% R = handles.Pipeline.OrigRed;
% G = handles.Pipeline.OrigGreen;
% B = handles.Pipeline.OrigBlue;
%I = uint16(cat(3,R,G,B));
%I = imread(strcat(path,imnames(imind,:),num2str(channel),'.tif'));
%L = imread([path imnames(imind,:) '0_SegmentedCells.tif']);
fieldname = ['Segmented',strObjectName];
L = uint16(CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale'));
G = handles.Pipeline.(strImageName);

%figure(imind)
%h(1) = subplot(221)
%imshow(L,[]);
%h(2) = subplot(222)
%imshow(log(double(I))/max(max(max(log(double(I))))),[]);

%%%%infotypeVAR01 = objectgroup

imind = handles.Current.SetBeingAnalyzed;
%channel = 1;
%stats = regionprops(L, 'PixelIdxList', 'PixelList');
%b = bwboundaries(L);
%centerpoints = handles.Measurements.Cells.Location{imind};

%figure(imind)
%h(3) = subplot(223)
%imshow(ves)
%linkaxes(h) 

% extract cell level features
[features,fnames] = extract_cell_features(double(L),G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

featurefieldname = [ModuleName '_',strImageName,'Features'];
fieldname = [ModuleName '_',strImageName];
%handles.Measurements.(objectname).(modulename_RescaledRed){imind}(n_objects, m_measurements)
%handles.Measurements.(objectname).(modulename_RescaledRedFeatures) = {'column1_description','column2_description'};
handles.Measurements.(strObjectName).(fieldname){imind} = features;
handles.Measurements.(strObjectName).(featurefieldname) = fnames;


end
