function plot_manopt_examples(results, family_filter)
%PLOT_MANOPT_EXAMPLES One PDF per Manopt-derived pilot, using seed 1.
% Optional family_filter selects which existing results to redraw.
startup_rhzls;
root=fileparts(mfilename('fullpath')); folder=fullfile(root,'docs');
if nargin<1
    data=load(fullfile(folder,'manopt_candidate_results.mat'),'results'); results=data.results;
end
families={'generalized_eigen','thomson_tetrahedron','maxcut_cycle5'};
if nargin < 2, family_filter = families; end
if ischar(family_filter) || isstring(family_filter), family_filter = cellstr(family_filter); end
assert(iscellstr(family_filter) && all(ismember(family_filter,families)), ...
    'Select one or more supported Manopt-derived families.');
titles={'Generalized eigenproblem / Grassmann','Thomson tetrahedron / product of spheres', ...
    'Max-Cut SDP / cycle $C_5$'};
for k=1:numel(families)
    if ~ismember(families{k},family_filter), continue; end
    entry=results(find(strcmp({results.family},families{k}) & [results.seed]==1,1));
    if ~isempty(entry.manopt.error) || ~isempty(entry.rhzls.error)
        warning('rhzls:CandidateFailed','Skipping failed candidate %s.',families{k}); continue;
    end
    label='$\|X_k^TX_k-G_\star\|_F/\|G_\star\|_F$';
    if strcmp(families{k},'generalized_eigen')
        label='$\|P_k-P_\star\|_F/\|P_\star\|_F$';
    end
    comparison=struct('name',families{k},'display_title',titles{k}, ...
        'rhzls_label','RHZLS (approx. Wolfe)','manopt',entry.manopt,'rhzls',entry.rhzls);
    comparison.reference=struct('optimum',entry.optimum,'solution_label',label);
    [fig,curves]=rhzls.plot_comparison(comparison);
    assert(numel(curves)==6);
    assert(~any(contains(fig.CurrentAxes.Legend.String,'either')));
    rhzls.export_pdf(fig,fullfile(folder,['candidate_' families{k} '.pdf']));
    if strcmp(char(fig.Visible),'off'), close(fig); end
end
end
