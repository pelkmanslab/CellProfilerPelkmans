function [] = checkTestInit()
% checks that environment variables and initial working directory are
% correct

    if isempty(getenv('CELLPROFILER_PELKMANS_TESTS_DIR'))
        error('Please set the environmental variable CELLPROFILER_PELKMANS_TESTS_DIR with the location of test data');
    end


    % Check that we are in right directory
    [upperPath, deepestFolder, ignoreThisStr] = fileparts(pwd());
    if (strcmp(deepestFolder,'Scripts')==0)
       error('Not in Scripts/ directory');
    end
end
