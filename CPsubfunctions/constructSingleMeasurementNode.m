function measurementNode = constructSingleMeasurementNode(measurementNode,measurementData, measurementName, measuerementFunction,fieldname)

    if isfield(measurementNode,'measurementField')<1
        measurementNode.measurementField = cell(0);
    end
    
    if isfield(measurementNode,'measurementName')<1
        measurementNode.(measurementName).fieldNames = cell(0);
        measurementNode.measurementField{end+1} = measurementName;
    else
       if isfield( measurementNode.(measurementName), fieldNames )<1
           measurementNode.(measurementName).fieldNames = cell(0); 
       end
    end
    
    measurementNode.(measurementName).fieldNames = cell(0);
    measurementNode.(measurementName).( [fieldname 'Value'] ) = feval(measuerementFunction, measurementData);
    measurementNode.(measurementName).fieldNames{end+1} = fieldname;
