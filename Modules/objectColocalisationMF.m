function handles = objectColocalisationMF(handles)

% Help for the objectColocalisation module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% MatF in Pelkmans Lab, 010911.
% this module is analysing the colocalisation of object 3 with object 1 and  with object 2.
% these objects must be specified in the GUI
% Beta version 1


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Object 1 please?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
Item1Name = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Object 2 please?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
Item2Name = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Object 3 please?
%infotypeVAR03 = objectgroup
%inputtypeVAR03 = popupmenu
Item3Name = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What object is the parent?
%infotypeVAR04 = objectgroup
%inputtypeVAR04 = popupmenu
ParentName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Retrieves the label matrix image that contains the 
%%% segmented objects 1.
img_Item1 = CPretrieveimage(handles,['Segmented', Item1Name],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the 
%%% segmented objects 2.
img_Item2 = CPretrieveimage(handles,['Segmented', Item2Name],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the 
%%% segmented objects 3.
img_Item3 = CPretrieveimage(handles,['Segmented', Item3Name],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the parents
%%% segmented objects.
img_Parents = CPretrieveimage(handles,['Segmented', ParentName],ModuleName,'MustBeGray','DontCheckScale');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Colocalisation calculation%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Logicalimg_Item1=(img_Item1 >= 1);
Logicalimg_Item2=(img_Item2 >= 1);
Logicalimg_Item3=(img_Item3 >= 1);

max_ParentLabel = max(max(img_Parents));

%%% let's give to Objects  the label of their parents

tmpMatrix = double(Logicalimg_Item1) ;
tmpMatrix(tmpMatrix==1)= max_ParentLabel+1 ;
mask = (img_Parents - tmpMatrix) > 0 ;

img_Item1LabeledPerCell = img_Parents ;
img_Item1LabeledPerCell(mask) = 0 ; 

tmpMatrix = double(Logicalimg_Item2) ;
tmpMatrix(tmpMatrix==1)= max_ParentLabel+1 ;
mask = (img_Parents - tmpMatrix) > 0 ;

img_Item2LabeledPerCell = img_Parents ;
img_Item2LabeledPerCell(mask) = 0 ;

tmpMatrix = double(Logicalimg_Item3) ;
tmpMatrix(tmpMatrix==1)= max_ParentLabel+1 ;
mask = (img_Parents - tmpMatrix) > 0 ;

img_Item3LabeledPerCell = img_Parents ;
img_Item3LabeledPerCell(mask) = 0 ;

%%%let's calculate the surface of the Objects in each parents 

matSurfaceOfItem1PerCell = regionprops(img_Item1LabeledPerCell,'area');
matSurfaceOfItem1PerCell = cat(1,matSurfaceOfItem1PerCell.Area);

%you'll find this loop all over the code, it is a security check as missing labels
%at the end of the list are not counted by region props.
%
%ex: if you regionprops objects that are only present in cells 1 2 3 8 10
%(and thus carrying the same label), automatically 4 5 6 7 and 9
%get a 0, but let say, in that case, there are 12 cells, 11 and 12 are w/o
%objects;
%these two objects don't get a zero, they are just not taken into account,
%this is of course the mess for ratio-ing over arrays (we need to keep always all the labels). So I fix
%this with this loop.

if length(matSurfaceOfItem1PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfItem1PerCell)+1):max_ParentLabel);
    matSurfaceOfItem1PerCell(missinglabels)= 0;
end
    
matSurfaceOfItem2PerCell = regionprops(img_Item2LabeledPerCell,'area');
matSurfaceOfItem2PerCell = cat(1,matSurfaceOfItem2PerCell.Area);

if length(matSurfaceOfItem2PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfItem2PerCell)+1):max_ParentLabel);
    matSurfaceOfItem2PerCell(missinglabels)= 0;
end

matSurfaceOfItem3PerCell = regionprops(img_Item3LabeledPerCell,'area');
matSurfaceOfItem3PerCell = cat(1,matSurfaceOfItem3PerCell.Area);

if length(matSurfaceOfItem3PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfItem3PerCell)+1):max_ParentLabel);
    matSurfaceOfItem3PerCell(missinglabels)= 0;
end

%%%let's calculate the overlap between Objects
%%%Object1 and 3

Logicalimg_ColocItem1And3 = double((Logicalimg_Item1+Logicalimg_Item3)==2);
Logicalimg_ColocItem1And3(Logicalimg_ColocItem1And3==1)= max_ParentLabel+1 ;
mask = (img_Parents - Logicalimg_ColocItem1And3) > 0 ;

img_ColocItem1And3_LabeledPerCell = img_Parents ;
img_ColocItem1And3_LabeledPerCell(mask) = 0 ;

matSurfaceOfColocItem1And3PerCell = regionprops(img_ColocItem1And3_LabeledPerCell,'area');
matSurfaceOfColocItem1And3PerCell = cat(1,matSurfaceOfColocItem1And3PerCell.Area);

if length(matSurfaceOfColocItem1And3PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfColocItem1And3PerCell)+1):max_ParentLabel);
    matSurfaceOfColocItem1And3PerCell(missinglabels)= 0;
end

%%%Object 2 and 3

Logicalimg_ColocItem2And3 = double((Logicalimg_Item2+Logicalimg_Item3)==2);
Logicalimg_ColocItem2And3(Logicalimg_ColocItem2And3==1)= max_ParentLabel+1 ;
mask = (img_Parents - Logicalimg_ColocItem2And3) > 0 ;

img_ColocItem2And3_LabeledPerCell = img_Parents ;
img_ColocItem2And3_LabeledPerCell(mask) = 0 ;

matSurfaceOfColocItem2And3PerCell = regionprops(img_ColocItem2And3_LabeledPerCell,'area');
matSurfaceOfColocItem2And3PerCell = cat(1,matSurfaceOfColocItem2And3PerCell.Area);

if length(matSurfaceOfColocItem2And3PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfColocItem2And3PerCell)+1):max_ParentLabel);
    matSurfaceOfColocItem2And3PerCell(missinglabels)= 0;
end

%%%Object 1 and 2 as well as common part of 1 and 2 with 3

Logicalimg_ColocItem1And2 = double((Logicalimg_Item1+Logicalimg_Item2)==2);
mask4overlapping1and2 = logical(Logicalimg_ColocItem1And2);

Logicalimg_ColocItem12and3 = double((Logicalimg_ColocItem1And2+Logicalimg_Item3)==2);

Logicalimg_ColocItem1And2(Logicalimg_ColocItem1And2==1)= max_ParentLabel+1 ;
mask = (img_Parents - Logicalimg_ColocItem1And2) > 0 ;

img_ColocItem1And2_LabeledPerCell = img_Parents ;
img_ColocItem1And2_LabeledPerCell(mask) = 0 ;

matSurfaceOfColocItem1And2PerCell = regionprops(img_ColocItem1And2_LabeledPerCell,'area');
matSurfaceOfColocItem1And2PerCell = cat(1,matSurfaceOfColocItem1And2PerCell.Area);

if length(matSurfaceOfColocItem1And2PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfColocItem1And2PerCell)+1):max_ParentLabel);
    matSurfaceOfColocItem1And2PerCell(missinglabels)= 0;
end

Logicalimg_ColocItem12and3(Logicalimg_ColocItem12and3==1)= max_ParentLabel+1 ;
mask = (img_Parents - Logicalimg_ColocItem12and3) > 0 ;

img_ColocItem12And3_LabeledPerCell = img_Parents ;
img_ColocItem12And3_LabeledPerCell(mask) = 0 ;

matSurfaceOfColocItem12And3PerCell = regionprops(img_ColocItem12And3_LabeledPerCell,'area');
matSurfaceOfColocItem12And3PerCell = cat(1,matSurfaceOfColocItem12And3PerCell.Area);

if length(matSurfaceOfColocItem12And3PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfColocItem12And3PerCell)+1):max_ParentLabel);
    matSurfaceOfColocItem12And3PerCell(missinglabels)= 0;
end


%%%The overlap between object 1 and 2 allows us to keep only the non
%%%overlapping lyso and golgi, let's do it now.

%calculation of the area of Object1 not colocalised with Object2 for each
%parent
nonOverlappingObject1 = img_Item1LabeledPerCell;
nonOverlappingObject1(mask4overlapping1and2) = 0 ; %keep only the non overlapping object 1, labelled per parent object
matSurfaceOfNonoverlappingItem1PerCell = regionprops(nonOverlappingObject1,'area');
matSurfaceOfNonoverlappingItem1PerCell = cat(1,matSurfaceOfNonoverlappingItem1PerCell.Area);

if length(matSurfaceOfNonoverlappingItem1PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfNonoverlappingItem1PerCell)+1):max_ParentLabel);
    matSurfaceOfNonoverlappingItem1PerCell(missinglabels)= 0;
end

%calculation of the area of Object2 not colocalised with Object1 for each
%parent
nonOverlappingObject2 = img_Item2LabeledPerCell;
nonOverlappingObject2(mask4overlapping1and2) = 0 ; %keep only the non overlapping object 2, labelled per parent object
matSurfaceOfNonoverlappingItem2PerCell = regionprops(nonOverlappingObject2,'area');
matSurfaceOfNonoverlappingItem2PerCell = cat(1,matSurfaceOfNonoverlappingItem2PerCell.Area);

if length(matSurfaceOfNonoverlappingItem2PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfNonoverlappingItem2PerCell)+1):max_ParentLabel);
    matSurfaceOfNonoverlappingItem2PerCell(missinglabels)= 0;
end

%calculation of the area of non overlapping Object1 colocalising with Object3 for each
%parent
logical_nonOverlappingObject1 = logical(nonOverlappingObject1);
img_colocNonOverlappingObject1with3 = double((logical_nonOverlappingObject1+Logicalimg_Item3)==2);
img_colocNonOverlappingObject1with3(img_colocNonOverlappingObject1with3 == 1) = max_ParentLabel+1; 
mask = (img_Parents - img_colocNonOverlappingObject1with3) > 0 ;
img_colocNonOverlappingObject1with3_LabeledPerCell = img_Parents ;
img_colocNonOverlappingObject1with3_LabeledPerCell(mask) = 0 ;

matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell = regionprops(img_colocNonOverlappingObject1with3_LabeledPerCell,'area');
matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell = cat(1,matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell.Area);

if length(matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell)+1):max_ParentLabel);
    matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell(missinglabels)= 0;
end


%calculation of the area of non overlapping Object2 colocalising with Object3 for each
%parent
logical_nonOverlappingObject2 = logical(nonOverlappingObject2);
img_colocNonOverlappingObject2with3 = double((logical_nonOverlappingObject2+Logicalimg_Item3)==2);
img_colocNonOverlappingObject2with3(img_colocNonOverlappingObject2with3 == 1) = max_ParentLabel+1; 
mask = (img_Parents - img_colocNonOverlappingObject2with3) > 0 ;
img_colocNonOverlappingObject2with3_LabeledPerCell = img_Parents ;
img_colocNonOverlappingObject2with3_LabeledPerCell(mask) = 0 ;

matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell = regionprops(img_colocNonOverlappingObject2with3_LabeledPerCell,'area');
matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell = cat(1,matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell.Area);

if length(matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell)<max_ParentLabel;
    missinglabels  = ((length(matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell)+1):max_ParentLabel);
    matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell(missinglabels)= 0;
end


%%%Let's do some ratio
%%%
%%%fraction of Object 1 occupied by Object 3
mat_FractionOfObject1_OccupiedPerObject3 = matSurfaceOfColocItem1And3PerCell ./ matSurfaceOfItem1PerCell ;
%%%fractino of Object 2 occupied by Object 3
mat_FractionOfObject2_OccupiedPerObject3 = matSurfaceOfColocItem2And3PerCell ./ matSurfaceOfItem2PerCell ;
%%%fraction of Object 3 occupied by Object 1
mat_FractionOfObject3_OccupiedPerObject1 = matSurfaceOfColocItem1And3PerCell ./ matSurfaceOfItem3PerCell ;
%%%fraction of Object 3 occupied by Object 2
mat_FractionOfObject3_OccupiedPerObject2 = matSurfaceOfColocItem2And3PerCell ./ matSurfaceOfItem3PerCell ;
%%%fraction of object 1 non overlapping with object 2
mat_FractionOfremainingObj1 = matSurfaceOfNonoverlappingItem1PerCell ./ matSurfaceOfItem1PerCell;
%%%fraction of object 2 non overlapping with object 1
mat_FractionOfremainingObj2 = matSurfaceOfNonoverlappingItem2PerCell ./ matSurfaceOfItem2PerCell;
%%%fraction of non overlapping Object 1 occupied by Object 3
mat_FractionOfremainingObj1__OccupiedPerObject3 = matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell ./ matSurfaceOfNonoverlappingItem1PerCell;
%%%fraction of non overlapping Object 2 occupied by Object 3
mat_FractionOfremainingObj2__OccupiedPerObject3 = matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell ./ matSurfaceOfNonoverlappingItem2PerCell;

mat_SurfaceOfObj3inALLNO = matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell + matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell;

mat_FractionOfremainingObj3inAllNO_inObj1 = matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell/mat_SurfaceOfObj3inALLNO;
mat_FractionOfremainingObj3inAllNO_inObj2 = matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell/mat_SurfaceOfObj3inALLNO;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%as a feature we have: Objects surface in each parent, Surface of
%%%collocalisation between objects in each parent and ratios 

handles.Measurements.Cells.objectsCollocalisationFeatures{1} = [Item1Name '_surface_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{2} = [Item2Name '_surface_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{3} = [Item3Name '_surface_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{4} = [Item3Name '_in_' Item1Name '_surface_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{5} = [Item3Name '_in_' Item2Name '_surface_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{6} = [Item1Name '_in_' Item2Name '_surface_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{7} = ['fraction_of_' Item1Name '_occupied_by_' Item3Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{8} = ['fraction_of_' Item2Name '_occupied_by_' Item3Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{9} = ['fraction_of_' Item3Name '_occupied_by_' Item1Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{10} = ['fraction_of_' Item3Name '_occupied_by_' Item2Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{11} = ['fraction_of_non_overlapping_' Item1Name '_with_' Item2Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{12} = ['fraction_of_non_overlapping_' Item2Name '_with_' Item1Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{13} = ['fraction_of_non_overlapping_' Item1Name '_occupied_by_' Item3Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{14} = ['fraction_of_non_overlapping_' Item2Name '_occupied_by_' Item3Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{15} = ['fraction_of_' Item3Name '_in_NOA_located_in_' Item1Name '_in_' ParentName];
handles.Measurements.Cells.objectsCollocalisationFeatures{16} = ['fraction_of_' Item3Name '_in_NOA_located_in_' Item2Name '_in_' ParentName];

matToStore = cat(2,matSurfaceOfItem1PerCell,matSurfaceOfItem2PerCell,matSurfaceOfItem3PerCell,...
    matSurfaceOfColocItem1And3PerCell,matSurfaceOfColocItem2And3PerCell,matSurfaceOfColocItem1And2PerCell,...
    mat_FractionOfObject1_OccupiedPerObject3, mat_FractionOfObject2_OccupiedPerObject3,mat_FractionOfObject3_OccupiedPerObject1,mat_FractionOfObject3_OccupiedPerObject2,...
    mat_FractionOfremainingObj1, mat_FractionOfremainingObj2, mat_FractionOfremainingObj1__OccupiedPerObject3, mat_FractionOfremainingObj2__OccupiedPerObject3,...
    mat_FractionOfremainingObj3inAllNO_inObj1, mat_FractionOfremainingObj3inAllNO_inObj2);
    
handles.Measurements.Cells.objectsCollocalisation{handles.Current.SetBeingAnalyzed} = matToStore ;



clear ('matSurfaceOfItem1PerCell','matSurfaceOfItem2PerCell','matSurfaceOfItem3PerCell',...
    'matSurfaceOfColocItem1And3PerCell','matSurfaceOfColocItem2And3PerCell','matSurfaceOfColocItem1And2PerCell',...
    'mat_FractionOfObject1_OccupiedPerObject3', 'mat_FractionOfObject2_OccupiedPerObject3','mat_FractionOfObject3_OccupiedPerObject1','mat_FractionOfObject3_OccupiedPerObject2',...
    'mat_FractionOfremainingObj1', 'mat_FractionOfremainingObj2', 'mat_FractionOfremainingObj1__OccupiedPerObject3', 'mat_FractionOfremainingObj2__OccupiedPerObject3',...
    'matSurfaceOfNonoverlappingItem1PerCell','matSurfaceOfNonoverlappingItem2PerCell','matSurfaceOfNonoverlappingItem1_overlappingWithObject3_PerCell',...
    'matSurfaceOfNonoverlappingItem2_overlappingWithObject3_PerCell', 'matSurfaceOfColocItem12And3PerCell','matSurfaceOfColocItem1And2PerCell' );
end

