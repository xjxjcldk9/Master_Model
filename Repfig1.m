%%%%%%Replicate fig1
global oo_ M_ options_ % get Dynare structures;


n = 60;
output = zeros(length(M_.endo_names), n);

options_.noprint=1; %關掉通知
options_.order=2;
options_.irf=0;

set_param_value('Me', 1e-10);
set_param_value('Phib', 0.0);



LTV = linspace(1e-10,0.96,n);
for i = 1:n
    set_param_value('Me', LTV(i));
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, M_.endo_names);
    output(:,i) = oo_.mean(:,1);
end

save('Fig1original', 'output');

set_param_value('Alpha', 0.1);
for i = 1:n
    set_param_value('Me', LTV(i));
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, M_.endo_names);
    output(:,i) = oo_.mean(:,1);
end
save('Fig1alpha', 'output');

set_param_value('Alpha', 0.7926);
