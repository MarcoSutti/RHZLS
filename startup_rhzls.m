function paths = startup_rhzls(manopt_root)
%STARTUP_RHZLS Add this project and one Manopt source tree to the path.
%   startup_rhzls() uses the bundled Manopt, or the last explicit selection.
%   startup_rhzls('C:\path\to\manopt') selects a different Manopt checkout.
%   The argument may name the checkout root or its inner manopt directory.
%   Only this project root and the selected Manopt sources are added.

    persistent selected_root selected_source selected_version_root
    project_root = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(manopt_root)
        if isempty(selected_root)
            manopt_root = fullfile(project_root, 'external', 'manopt');
        else
            manopt_root = selected_root;
        end
    end
    if ~(ischar(manopt_root) || (isstring(manopt_root) && isscalar(manopt_root)))
        error('rhzls:InvalidManoptRoot', 'Manopt root must be a folder path.');
    end
    manopt_root = char(manopt_root);
    if exist(fullfile(manopt_root, 'manopt', 'core'), 'dir') == 7
        source_root = fullfile(manopt_root, 'manopt');
    elseif exist(fullfile(manopt_root, 'core'), 'dir') == 7
        source_root = manopt_root;
    else
        error('rhzls:MissingManopt', ...
            ['Manopt sources not found in "%s". Restore external/manopt, ' ...
             'or call startup_rhzls with your Manopt checkout path.'], manopt_root);
    end
    [~, root_info] = fileattrib(manopt_root);
    [~, source_info] = fileattrib(source_root);
    manopt_root = root_info.Name;
    source_root = source_info.Name;
    checkout_root = manopt_root;
    if strcmpi(manopt_root, source_root) ...
            && exist(fullfile(fileparts(source_root), 'manopt_version.m'), 'file') == 2
        checkout_root = fileparts(source_root);
    end

    if ~isempty(selected_source) && ~strcmpi(selected_source, source_root)
        rmpath(genpath(selected_source));
    end
    if ~isempty(selected_version_root) && ~strcmpi(selected_version_root, checkout_root)
        rmpath(selected_version_root);
        selected_version_root = [];
    end
    addpath(genpath(source_root), '-begin');
    % The official checkout keeps version and license files at its root.
    if exist(fullfile(checkout_root, 'manopt_version.m'), 'file') == 2
        addpath(checkout_root, '-begin');
        selected_version_root = checkout_root;
    end
    addpath(project_root, '-begin');
    required = {'spherefactory', 'stiefelfactory', 'fixedrankembeddedfactory', ...
        'steepestdescent', 'linesearch'};
    for k = 1:numel(required)
        if exist(required{k}, 'file') ~= 2
            error('rhzls:MissingManopt', 'Manopt function %s is missing.', required{k});
        end
        resolved_file = which(required{k});
        if ~startsWith(lower(resolved_file), lower([source_root filesep]))
            error('rhzls:ManoptPathConflict', ...
                'Function %s resolves outside selected Manopt: %s', required{k}, resolved_file);
        end
    end
    selected_root = manopt_root;
    selected_source = source_root;
    paths = struct('project', project_root, 'manopt', source_root);
    if exist(fullfile(checkout_root, 'manopt_version.m'), 'file') == 2
        [paths.manopt_version, paths.manopt_release] = manopt_version();
    end
    manifest_file = fullfile(checkout_root, 'UPSTREAM_MANIFEST.json');
    if exist(manifest_file, 'file') == 2
        provenance = jsondecode(fileread(manifest_file));
        paths.manopt_commit = provenance.commit;
    end
end
