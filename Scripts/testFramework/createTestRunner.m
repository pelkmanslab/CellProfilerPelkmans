function [ runner ] = createTestRunner( outFilePathRelative )
% creates a TestRunner that writes a TAP file to outFilePathRelative

    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.TAPPlugin
    import matlab.unittest.plugins.ToFile

    runner = TestRunner.withTextOutput;

    % NOTE: It is critically important that an absolute file path is passed to
    %   TAPPlugin.producingOriginalFormat as otherwise test entries become
    %   missing in the output TAP file
    outFileResolved = fullfile(pwd(),outFilePathRelative);

    % We delete the existing tapFile as TAPPlugin.producingOriginalFormat
    %   appends to any existing files, and we need it to be only a s
    %    single test case for Jenkins.
    delete(outFileResolved);
    plugin = TAPPlugin.producingOriginalFormat(ToFile(outFileResolved));

    runner.addPlugin(plugin)

end