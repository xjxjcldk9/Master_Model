
%%%%%%Fig2
global oo_ M_ options_ % get Dynare structures;


options_.noprint=1; %關掉通知
options_.order=2;
options_.irf=0;


set_param_value('Me', 0.90);
set_param_value('Phib', 0.0);

set_param_value('Rpi', 1.1955);
set_param_value('Rs', -0.0968);
set_param_value('Ry', 0.1707);

set_param_value('Rq', 0.0);
set_param_value('Phiq', 0.0);
set_param_value('Phiy', 0.0);
set_param_value('Phib', 0.0);


[info, oo_, options_] = stoch_simul(M_, options_, oo_, {'Ws', 'Wb', 'We'});


wws = oo_.mean(1);
wwb = oo_.mean(2);
wwe = oo_.mean(3);

wtotp = (1-Betas)*wws + (1-Betab)*wwb + (1-Betae)*wwe;

n=100;
Phibb = linspace(0,3.5,n);
for i = 1:n
    set_param_value('Phib', Phibb(i));

    [info, oo_, options_] = stoch_simul(M_, options_, oo_, {'Ws', 'Wb', 'We'});
    Wsmp = oo_.mean(1);
    Wbmp = oo_.mean(2);
    Wemp = oo_.mean(3);

    CEs(i) = exp((1-Betas)*(Wsmp - wws)) - 1;
    CEb(i) = exp((1-Betab)*(Wbmp - wwb)) - 1;
    CEe(i) = exp((1-Betae)*(Wemp - wwe)) - 1;

    wtotap = (1-Betas)*Wsmp + (1-Betab)*Wbmp + (1-Betae)*Wemp;
    CE(i) = exp(1/3*(wtotap - wtotp)) - 1;
end

plot(Phibb,CEs,'--',Phibb,CEb,'--',Phibb,CEe,':',Phibb,CE,'LineWidth',3);
legend('CEs','CEb','CEe','CE');

