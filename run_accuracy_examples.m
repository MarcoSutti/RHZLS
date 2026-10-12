function results = run_accuracy_examples(draw_plot)
%RUN_ACCURACY_EXAMPLES Fixed benchmark suite; report all seeds and failures.
% Pass false to save results without drawing or exporting intermediate PDFs.
if nargin < 1, draw_plot = true; end
validateattributes(draw_plot, {'logical'}, {'scalar'});
manopt_source = startup_rhzls;
root = fileparts(mfilename('fullpath'));
state = rng; cleanup = onCleanup(@() rng(state)); %#ok<NASGU>
families = {'sphere_positive','stiefel_positive','orthogonal_procrustes','rank_approximation'};
results = struct([]);
for family = families
    for seed = 1:3
        rng(seed,'twister');
        [problem,x0,fstar,target] = build_example(family{1});
        options = struct('maxiter',1000,'tolgradnorm',0,'minstepsize',0,'verbosity',0);
        if strcmp(family{1}, 'sphere_positive')
            options.maxiter = 200;
        elseif strcmp(family{1}, 'stiefel_positive')
            options.maxiter = 850;
        elseif strcmp(family{1}, 'orthogonal_procrustes')
            options.maxiter = 60;
        elseif strcmp(family{1}, 'rank_approximation')
            options.maxiter = 250;
        end
        options.statsfun = @(p,x,s) record_stats(p,x,s,fstar,target);
        entry = struct('family',family{1},'seed',seed,'problem',problem,'x0',x0, ...
            'optimum',fstar,'target',target);
        for method = {'manopt','rhzls','rhzls_approx'}
            name = method{1};
            method_options = options;
            if ~strcmp(name,'manopt')
                method_options.linesearch = @rhzls.manopt_linesearch;
                method_options.rhzls = struct('acceptance','either');
                if strcmp(name,'rhzls_approx')
                    method_options.rhzls.acceptance = 'approximate';
                end
            end
            try
                [x,cost,info,effective_options] = steepestdescent(problem,x0,method_options);
                if strcmp(name,'manopt')
                    assert(isequal(effective_options.linesearch,@linesearch));
                end
                info = rhzls.add_evaluation_counts(info,effective_options.linesearch);
                entry.(name) = struct('x',x,'cost',cost,'info',info,'options',effective_options,'error','');
                fprintf('%s seed=%d %s: iter=%d grad=%.12e solution=%.12e\n', ...
                    family{1},seed,name,info(end).iter,info(end).gradnorm,info(end).solution_error);
            catch exception
                entry.(name) = struct('x',[],'cost',NaN,'info',[],'options',method_options,'error',exception.message);
                fprintf('%s seed=%d %s FAILED: %s\n',family{1},seed,name,exception.message);
            end
        end
        if isempty(results), results = entry;
        else, results(end+1) = entry; end %#ok<AGROW>
    end
end
save(fullfile(root,'docs','accuracy_examples.mat'),'results','manopt_source');
if draw_plot, plot_accuracy_examples(results); end
end

function [p,x0,fstar,target] = build_example(family)
switch family
    case 'sphere_positive'
        A = diag(linspace(1,5,12));
        [p,x0,ref] = rhzls.make_problem('sphere',struct('A',A));
        fstar = ref.optimum; target = eye(12,1);
    case 'stiefel_positive'
        A = diag(linspace(1,5,12));
        [p,x0,ref] = rhzls.make_problem('stiefel',struct('A',A,'p',2,'weights',[2 1]));
        fstar = ref.optimum; target = eye(12,2);
    case 'orthogonal_procrustes'
        [M,curve,maxstep] = rhzls.manifold('stiefel',[4 4]);
        A = diag(1:4);
        p = struct('M',M,'curve',curve,'maxstep',maxstep);
        p.cost = @(x) .5*norm(x-A,'fro')^2;
        p.egrad = @(x) x-A;
        [x0,~] = qr(randn(4));
        if det(x0)<0, x0(:,1)=-x0(:,1); end
        target = eye(4); fstar = .5*norm(target-A,'fro')^2;
    case 'rank_approximation'
        [M,curve,maxstep] = rhzls.manifold('fixedrank',[12 12 2]);
        spectrum = linspace(2,.2,12); A = diag(spectrum);
        p = struct('M',M,'curve',curve,'maxstep',maxstep);
        p.cost = @(x) .5*norm(rhzls.full_matrix(x)-A,'fro')^2;
        p.egrad = @(x) rhzls.full_matrix(x)-A;
        [U,~] = qr(eye(12,2)+.1*randn(12,2),0);
        [V,~] = qr(eye(12,2)+.1*randn(12,2),0);
        x0 = struct('U',U,'S',.8*diag(spectrum(1:2)),'V',V);
        target = diag([spectrum(1:2),zeros(1,10)]);
        fstar = .5*sum(spectrum(3:end).^2);
end
end

function stats = record_stats(problem,x,stats,fstar,target)
stats.objective_error = abs(stats.cost-fstar);
X = rhzls.full_matrix(x);
if isvector(target)
    stats.solution_error = min(norm(X-target),norm(X+target));
elseif size(target,1) ~= size(target,2)
    signs = sign(sum(target.*X,1)); signs(signs==0)=1;
    stats.solution_error = norm(X-target.*signs,'fro')/norm(target,'fro');
else
    stats.solution_error = norm(X-target,'fro')/norm(target,'fro');
end
end
