function validate_plateau_iterations()
%VALIDATE_PLATEAU_ITERATIONS Check actual saved histories, without padding.
root=fileparts(mfilename('fullpath'));
a=load(fullfile(root,'accuracy_examples.mat'),'results');
families={'sphere_positive','stiefel_positive','orthogonal_procrustes','rank_approximation'};
limits=[200 850 60 250]; count=0;
for k=1:numel(a.results)
    budget=limits(strcmp(families,a.results(k).family));
    assert(isequal(a.results(k).manopt.options.linesearch,@linesearch));
    assert(all([a.results(k).manopt.info.ng]==0));
    for m={'manopt','rhzls','rhzls_approx'}
        h=a.results(k).(m{1}).info;
        assert(isequal([h.iter],0:budget)); count=count+1;
    end
end
c=load(fullfile(root,'manopt_candidate_results.mat'),'results');
for k=1:numel(c.results)
    assert(isequal(c.results(k).manopt.options.linesearch,@linesearch));
    assert(all([c.results(k).manopt.info.ng]==0));
    budget=1000;
    if strcmp(c.results(k).family,'maxcut_cycle5'), budget=60;
    elseif strcmp(c.results(k).family,'thomson_tetrahedron'), budget=50;
    end
    for m={'manopt','rhzls'}
        h=c.results(k).(m{1}).info;
        assert(isequal([h.iter],0:budget)); count=count+1;
    end
end
d=load(fullfile(root,'demo_results.mat'),'details');
for family={'sphere','stiefel'}
    c=d.details.(family{1}).comparison;
    assert(c.common_options.tolgradnorm==0);
    assert(isequal(c.manopt.details.options.linesearch,@linesearch));
    assert(all([c.manopt.info.ng]==0));
    for m={'manopt','rhzls'}
        h=c.(m{1}).info;
        assert(isequal([h.iter],0:c.common_options.maxiter)); count=count+1;
    end
end
assert(count==58);
fprintf('PASS: all 58 saved runs contain every iteration through maxiter.\n');
fprintf('PASS: every baseline run uses official @linesearch with ng=0.\n');
end
