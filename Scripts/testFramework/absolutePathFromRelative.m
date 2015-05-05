function [ absolutePath ] = absolutePathFromRelative( relativePath )
% generates an absolute path to the testing folder from a relative path

   basePath = getenv('CELLPROFILER_PELKMANS_TESTS_DIR');
   absolutePath = fullfile(basePath,relativePath);
end

