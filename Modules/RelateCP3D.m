function handles = RelateCP3D(handles)

% Help for the RelateCP3D module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Assigns relationships: All objects (e.g. speckles) within a parent object
% (e.g. nucleus) become its children.
% *************************************************************************
%
% Allows counting the number of objects within each parent object.
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What objects are the children objects (subobjects)?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
SubObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What are the parent objects?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
ParentName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which relate method do you want to apply?
%choiceVAR03 = Centroid
iMethod = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

fieldname = ['Segmented', SubObjectName];
if isfield(handles.Pipeline.(fieldname),'Format')
    switch handles.Pipeline.(fieldname).Format
        case 'SegmentationCC'
            SubObjectLabel = handles.Pipeline.(fieldname).Label;
        otherwise
            error(['Image processing was canceled in the ', ModuleName, ' module because Segmentation of ' SubObjectName ' follows undefined format.'])
    end
else
    SubObjectLabel =  handles.Pipeline.(fieldname);
end


fieldname = ['Segmented', ParentName];
if isfield(handles.Pipeline.(fieldname),'Format')
    switch handles.Pipeline.(fieldname).Format
        case 'SegmentationCC'
            ParentLabel = handles.Pipeline.(fieldname).Label;
        otherwise
            error(['Image processing was canceled in the ', ModuleName, ' module because Segmentation of ' SubObjectName ' follows undefined format.'])
    end
else
    ParentLabel =  handles.Pipeline.(fieldname);
end



%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

[numChildren ParentID]= relateobjectsCP3D(SubObjectLabel,ParentLabel,iMethod);

% Save Results
handles = CPaddmeasurements(handles,ParentName,'Children',[SubObjectName,'Count'],numChildren);
handles = CPaddmeasurements(handles,SubObjectName,'Parent',ParentName,ParentID);


end