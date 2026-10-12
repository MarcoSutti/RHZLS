function plot_accuracy_examples(results)
%PLOT_ACCURACY_EXAMPLES Normalized Armijo vs approximate-Wolfe RHZLS comparison.
startup_rhzls;
root = fileparts(mfilename('fullpath')); folder = fullfile(root,'docs');
if nargin < 1
    data = load(fullfile(folder,'accuracy_examples.mat'),'results'); results = data.results;
end
families = {'sphere_positive','stiefel_positive','orthogonal_procrustes','rank_approximation'};
titles = {'Positive Rayleigh / sphere','Positive Brockett / Stiefel', ...
    'Orthogonal Procrustes / SO(4)','Best rank-2 approximation'};
for k = 1:numel(families)
    entry = results(find(strcmp({results.family},families{k}) & [results.seed]==1,1));
    if ~isempty(entry.manopt.error) || ~isempty(entry.rhzls_approx.error), continue; end
    comparison = struct('name',families{k},'display_title',titles{k}, ...
        'manopt',entry.manopt,'rhzls',entry.rhzls_approx);
    comparison.rhzls_label = 'RHZLS (approx. Wolfe)';
    comparison.window_history = entry.rhzls_approx.info;
    comparison.reference = struct('optimum',entry.optimum, ...
        'solution_label','$d(x_k,\mathcal{X}_\star)/\|x_\star\|$');
    fig = rhzls.plot_comparison(comparison);
    rhzls.export_pdf(fig,fullfile(folder,['accuracy_' families{k} '.pdf']));
    if strcmp(char(fig.Visible),'off'), close(fig); end

end
end
