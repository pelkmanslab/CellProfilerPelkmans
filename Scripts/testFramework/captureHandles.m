function [handles] = captureHandles( handles )
%   Saves handles_in and handles_out before calling a delegate module
%   
%   To be used to create test-cases for modules
%    e.g. to create a test for DiscardObjectsBySize
%
%   In this case, the delgate is DiscardObjectsBySize. We inject
%    this function into normal execution.
%
%   Steps:
%   1. Temporarily rename the delegate. Update function name.
%      e.g. DiscardObjectsBySize.m -> DiscardObjectsBySizeTemp.m
%
%   2. Copy this function, and rename to the original delegate-name. Update
%      function name
%      e.g. captureHandles(copy) -> DiscardObjectsBySize.m
%
%   3. Rename the DelgateTemp call to DiscardObjectsBySizeTemp
%
%   4. Run CellProfiler, and each time the function is called, it will
%   write handles_in_001.mat etc. to the current working directory.
%
%   5. When finished, delete DiscardObjectsBySize.m and restore
%   DiscardObjectsBySizeTemp.m to its original name. Update function name.
%      DiscardObjectsBySizeTemp.m -> DiscardObjectsBySize
%
    persistent indexNum;
    
    if (isempty(indexNum))
       indexNum = 0;
    else
       indexNum = indexNum + 1;
    end

    save( sprintf('handles_in_%03d.mat',indexNum), handles );
    
    handles = DelegateTemp( handles );

    save( sprintf('handles_out_%03d.mat',indexNum), handles );
end

