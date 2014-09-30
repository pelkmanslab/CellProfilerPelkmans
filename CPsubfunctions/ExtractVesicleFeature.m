function handles = ExtractVesicleFeature(handles)

% Help for the ExtractVesicle module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Interface for vesicle detection by Sharif, Chowdhury
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4904 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object Name?
%defaultVAR01 = Vesicles
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Input Image Name In Matlab?
%defaultVAR02 = 
%infotypeVAR02 = objectgroup indep
inputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Cell Mask Name?
%defaultVAR03 = 
%infotypeVAR03 = objectgroup indep
maskName = char(handles.Settings.VariableValues{CurrentModuleNum,3});


%textVAR04 = Intial Cycle Counter?
%defaultVAR04 = 1
%infotypeVAR04 = objectgroup indep
cycleCounter = str2num( char( handles.Settings.VariableValues{CurrentModuleNum,4}) );


%textVAR05 = Riple's Starting Radius ?
%defaultVAR05 = 1
%infotypeVAR05 = objectgroup indep
rad1 = str2num( char( handles.Settings.VariableValues{CurrentModuleNum,5}) );

if rad1<1
    rad1 = 1;
end


%textVAR06 = Riple's Ending Radius (zero for disabling the feature) ?
%defaultVAR06 = 0
%infotypeVAR06 = objectgroup indep
rad2 = str2num( char( handles.Settings.VariableValues{CurrentModuleNum,6}) );




I = CPretrieveimage(handles,inputImageName,ModuleName,'MustBeGray','CheckScale');
fieldname = ['Segmented',maskName];
maskImage = CPretrieveimage(handles,fieldname,ModuleName,'DontCheckColor','DontCheckScale',size(I));
fieldname = ['Segmented', ObjectName];
vesicleImage = CPretrieveimage(handles,fieldname,ModuleName,'DontCheckColor','DontCheckScale',size(I));

[featmat,allvesicledata,fnames,vfnames, riplyMatrix] = extract_vesicle_features(maskImage,I,vesicleImage, rad1,rad2);


if isfield(handles.Measurements, 'cycleCount')<1
    handles.Measurements.('cycleCount') = cycleCounter; 
else
    handles.Measurements.cycleCount = handles.Measurements.cycleCount + 1;
end

handles.Measurements.('vesicleFeaturePerCell'){handles.Measurements.('cycleCount')}= featmat;
handles.Measurements.('IndividualVesicleData'){handles.Measurements.('cycleCount')}= allvesicledata;

handles.Measurements.('vesicleFeatureNamesPerCell' )  = fnames;
handles.Measurements.('IndividualVesicleFeatureNames') = vfnames;

if rad2 >=0 && rad1<=rad2
    handles.Measurements.('RiplesK'){handles.Measurements.('cycleCount')}= riplyMatrix;
end

end
