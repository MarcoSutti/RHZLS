function results = explore_manopt_examples(draw_plot)
%EXPLORE_MANOPT_EXAMPLES Reproducible pilots adapted from official Manopt examples.
% Pass false to save results without drawing or exporting intermediate PDFs.
% No existing benchmark or plotting output is replaced. Uses exact gradients.
% Sources:
% https://github.com/NicolasBoumal/manopt/blob/master/examples/generalized_eigenvalue_computation.m
% https://github.com/NicolasBoumal/manopt/blob/master/examples/thomson_problem.m
% https://github.com/NicolasBoumal/manopt/blob/master/examples/maxcut.m
if nargin < 1, draw_plot = true; end
validateattributes(draw_plot, {'logical'}, {'scalar'});
manopt_source = startup_rhzls;
root = fileparts(mfilename('fullpath'));
previous_rng = rng; restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
results = struct([]);
for family = {'generalized_eigen','thomson_tetrahedron','maxcut_cycle5'}
    for seed = 1:3
        rng(seed,'twister');
        [problem,x0,optimum,solutionfun] = build_problem(family{1});
        d = problem.M.randvec(x0);
        [~,v0] = problem.curve(x0,d,0);
        assert(norm(v0-d,'fro') < 1e-11);
        t = .15; h = 1e-6;
        [y,v] = problem.curve(x0,d,t);
        yp = problem.curve(x0,d,t+h); ym = problem.curve(x0,d,t-h);
        velocity_error = norm((yp-ym)/(2*h)-v,'fro')/max(1,norm(v,'fro'));
        egrad = problem.egrad(y);
        exact = sum(sum(egrad.*v));
        if isfield(problem, 'curve_slope')
            stable_slope = problem.curve_slope(y, egrad, v);
            assert(abs(stable_slope-exact) < 1e-10*max(1,abs(exact)));
        end
        finite_diff = (problem.cost(yp)-problem.cost(ym))/(2*h);
        slope_error = abs(exact-finite_diff)/max(1,abs(exact));
        assert(velocity_error < 1e-7 && slope_error < 1e-7);
        if isfield(problem, 'curve_slope')
            assert(abs(stable_slope-finite_diff)/max(1,abs(stable_slope)) < 1e-7);
        end
        if ~strcmp(family{1}, 'generalized_eigen')
            assert(max(abs(sum(y.*v,1))) < 1e-12*max(1,norm(v,'fro')));
            assert(max(abs(sum(y.^2,1)-1)) < 1e-12);
        end
        entry = struct('family',family{1},'seed',seed,'optimum',optimum, ...
            'x0',x0,'problem',problem,'velocity_error',velocity_error,'slope_error',slope_error);
        options = struct('maxiter',1000,'tolgradnorm',0,'minstepsize',0,'verbosity',0);
        if strcmp(family{1}, 'maxcut_cycle5')
            options.maxiter = 60;
        elseif strcmp(family{1}, 'thomson_tetrahedron')
            options.maxiter = 50;
        end
        options.statsfun = @(p,x,s) record_stats(x,s,optimum,solutionfun);
        for method = {'manopt','rhzls'}
            name = method{1};
            method_options = options;
            if ~strcmp(name,'manopt')
                method_options.linesearch = @rhzls.manopt_linesearch;
                method_options.rhzls = struct('acceptance','approximate');
            end
            try
                [x,cost,info,effective_options] = steepestdescent(problem,x0,method_options);
                if strcmp(name,'manopt')
                    assert(isequal(effective_options.linesearch,@linesearch));
                end
                info = rhzls.add_evaluation_counts(info,effective_options.linesearch);
                entry.(name) = struct('x',x,'cost',cost,'info',info,'options',effective_options,'error','');
                fprintf('%s seed=%d %s: iter=%d grad=%.12e relative solution=%.12e\n', ...
                    family{1},seed,name,info(end).iter,info(end).gradnorm,info(end).solution_error);
            catch exception
                entry.(name) = struct('x',[],'cost',NaN,'info',[],'options',method_options,'error',exception.message);
                fprintf('%s seed=%d %s FAILED: %s\n',family{1},seed,name,exception.message);
            end
        end
        if isempty(results), results=entry; else, results(end+1)=entry; end %#ok<AGROW>
    end
end
save(fullfile(root,'docs','manopt_candidate_results.mat'),'results','manopt_source');
if draw_plot, plot_manopt_examples(results); end
end

function [p,x0,optimum,solutionfun] = build_problem(family)
switch family
    case 'generalized_eigen'
        n=12; rank=2;
        B=diag(linspace(1,2,n)); L=chol(B,'lower');
        eigenvalues=linspace(1,5,n); A=L*diag(eigenvalues)*L';
        M=grassmanngeneralizedfactory(n,rank,B);
        curve=@(x,d,t) generalized_curve(x,d,t,L);
        p=struct('M',M,'curve',curve,'maxstep',@(x,d) Inf);
        p.M.retr=curve;
        p.cost=@(x) -sum(sum(x.*(A*x)));
        p.egrad=@(x) -2*A*x;
        % The exact curve velocity is horizontal. Pairing it with the
        % Riemannian gradient avoids cancellation by the normal component.
        p.curve_slope=@(x,g,v) M.inner(x,M.egrad2rgrad(x,g),v);
        [Q,~]=qr(randn(n,rank),0); x0=L'\Q;
        optimum=-sum(eigenvalues(end-rank+1:end));
        target=zeros(n); target(end-rank+1:end,end-rank+1:end)=eye(rank);
        solutionfun=@(x) norm((L'*x)*(L'*x)'-target,'fro')/sqrt(rank);
    case 'thomson_tetrahedron'
        n=4; M=obliquefactory(3,n);
        curve=@oblique_curve;
        p=struct('M',M,'curve',curve,'maxstep',@(x,d) Inf);
        p.M.retr=curve;
        p.cost=@thomson_cost;
        p.egrad=@thomson_gradient;
        % The exact normalization velocity is tangent to every sphere.
        % Remove the normal component before pairing to reduce cancellation.
        p.curve_slope=@(x,g,v) M.inner(x,M.egrad2rgrad(x,g),v);
        x0=M.rand(); optimum=9/32;
        target=(4/3)*eye(n)-ones(n)/3;
        solutionfun=@(x) norm(x'*x-target,'fro')/norm(target,'fro');
    case 'maxcut_cycle5'
        n=5; M=obliquefactory(2,n); curve=@oblique_curve;
        adjacency=zeros(n);
        for j=1:n, k=mod(j,n)+1; adjacency(j,k)=1; adjacency(k,j)=1; end
        laplacian=2*eye(n)-adjacency;
        p=struct('M',M,'curve',curve,'maxstep',@(x,d) Inf);
        p.M.retr=curve;
        p.cost=@(x) -.25*sum(sum(x.*(x*laplacian)));
        p.egrad=@(x) -.5*x*laplacian;
        p.curve_slope=@(x,g,v) M.inner(x,M.egrad2rgrad(x,g),v);
        x0=M.rand(); optimum=-2.5*(1+cos(pi/5));
        angles=4*pi*(0:n-1)/5;
        target=[cos(angles);sin(angles)]; target=target'*target;
        solutionfun=@(x) norm(x'*x-target,'fro')/norm(target,'fro');
end
end

function [y,v] = generalized_curve(x,d,t,L)
% Keep this entry point for function handles in previously saved MAT files.
if nargout>1
    [y,v]=rhzls.retractions.generalized_grassmann_curve(x,d,t,L);
else
    y=rhzls.retractions.generalized_grassmann_curve(x,d,t,L);
end
end

function [y,v] = oblique_curve(x,d,t)
% Keep this entry point for function handles in previously saved MAT files.
if nargout>1
    [y,v]=rhzls.retractions.oblique_curve(x,d,t);
else
    y=rhzls.retractions.oblique_curve(x,d,t);
end
end

function f = thomson_cost(x)
n=size(x,2); S=1-x'*x; indices=find(triu(true(n),1));
if any(S(indices)<=0), f=Inf; else, f=sum(1./S(indices))/n^2; end
end

function g = thomson_gradient(x)
n=size(x,2); S=1-x'*x; S(1:n+1:end)=Inf;
g=x*(1./S.^2)/n^2;
end

function s = record_stats(x,s,optimum,solutionfun)
s.objective_error=abs(s.cost-optimum);
s.solution_error=solutionfun(x);
end
