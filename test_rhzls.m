function report = test_rhzls(varargin)
%TEST_RHZLS Run scalar, retraction, pullback and optimization regressions.
%   TEST_RHZLS uses the Manopt located by STARTUP_RHZLS. An optional Manopt
%   directory can be supplied: TEST_RHZLS('path/to/manopt').
startup_rhzls(varargin{:});
test_directory = fullfile(fileparts(mfilename('fullpath')), 'tests');
old_path = path;
restore_path = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(test_directory);
report = run_rhzls_tests();
end
