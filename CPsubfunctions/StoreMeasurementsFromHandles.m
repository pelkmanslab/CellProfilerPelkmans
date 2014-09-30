function StoreMeasurementsFromHandles(handles)

% check if handles.Measurements exists
if not(isfield(handles, 'Measurements'))
   error('%s: Handles does not contain any Measurements',mfilename) 
end

strDataSorterOutputFolderPath = handles.Current.DefaultOutputDirectory;

%%%%%%%%%%%%%%%%%%%%%%%
%%%  STORING  DATA  %%%
%%%%%%%%%%%%%%%%%%%%%%%

ListOfObjects = fieldnames(handles.Measurements); 
for i = 1:length(ListOfObjects)
    ListOfMeasurements = fieldnames(handles.Measurements.(char(ListOfObjects(i))));
    for ii = 1:length(ListOfMeasurements)
        if isempty(regexp(ListOfMeasurements{ii},'.*Features$', 'once')) || isempty(regexp(ListOfMeasurements{ii},'.*Text$', 'once'))
            % we are not dealing with a ...Features or ...Text list

            OutPutFile = fullfile(strDataSorterOutputFolderPath, sprintf('Measurements_%s_%s',char(ListOfObjects(i)), char(ListOfMeasurements(ii))));
            matPossibleDescriptionIndex = [strcmp(ListOfMeasurements, strcat(char(ListOfMeasurements(ii)),'Features')) || strcmp(ListOfMeasurements, strcat(char(ListOfMeasurements(ii)),'Text'))];

            % initialize the current output variable Measurements
            Measurements = struct();            
            
            % add the current Measurement data
            Measurements.(char(ListOfObjects(i))).(char(ListOfMeasurements(intFieldDescriptionIndex))) = handles.Measurements.(char(ListOfObjects(i))).(char(ListOfMeasurements(intFieldDescriptionIndex)));
            
            % if this measurement field has a corresponding ...Features or
            % ...Text fieldname, include it in the Measurement structure.
            if ~isempty(find(matPossibleDescriptionIndex))
                Measurements.(char(ListOfObjects(i))).(char(ListOfMeasurements(ii))) = handles.Measurements.(char(ListOfObjects(i))).(char(ListOfMeasurements(ii)));
            end

            % store Measurement file
            disp(sprintf('%s: saving Measurements.%s.%s',mfilename,char(ListOfObjects(i)), char(ListOfMeasurements(ii))))            
            save(OutPutFile, 'Measurements');

            % clear Measurements for next round
            clear Measurements;
        else
            disp(sprintf('SKIPPED %s', char(ListOfMeasurements{ii})))
        end
    end
end


end % end function