function handles = CreateJobList(handles)

% Help for the Create Batch Files module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Produces YAML description of single jobs to run for CellProfilerPelkmans.
% *************************************************************************
%
% This module creates should be put at the end of any pipeline that 
% requires parallel (batched) image processing.
%
% Before using this module, you should read Help -> Getting Started ->
% Batch Processing. That help file also will instruct you on how to
% actually run the batch files that are created by this module.
%
% This code is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
% 
% Authors:
%     Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
% 
% Copyright 2014 Pelkmans group https://www.pelkmanslab.org
%
% $Revision: 1 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Which format is used to store job list?
%choiceVAR01 = YAML
FormatChoice = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu
%defaultVAR01 = YAML

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% If this isn't the first cycle, we are running on the
%%% cluster, and should just continue.
if (handles.Current.SetBeingAnalyzed > 1),
    return;
end

%%% Checks that this is the last module in the analysis path.
if (CurrentModuleNum ~= handles.Current.NumberOfModules),
    error(['Image processing was canceled because ', ModuleName, ' must be the last module in the pipeline.']);
end

%%% Do actual work.
cellpro.jobs.createList(handles);


CPhelpdlg('Job list been written. This analysis pipeline will now stop. You should proceed with cellpro --help executable to analyze jobs in parallel.', 'BatchFilesDialog');

%%% This is the first cycle, so this is the first time seeing this
%%% module.  It should cause a cancel so no further processing is done
%%% on this machine.
set(handles.timertexthandle,'string','Cancel')

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The figure window display is unnecessary for this module, so it is
%%% closed during the starting image cycle.
if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        close(ThisModuleFigureNumber)
    end
end
