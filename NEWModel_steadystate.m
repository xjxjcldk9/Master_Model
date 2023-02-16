function [ys,params,check] = NK_baseline_steadystate(ys,exo,M_,options_)
% function [ys,params,check] = NK_baseline_steadystate(ys,exo,M_,options_)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs:
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output:
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% Copyright (C) 2013-2020 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = M_.param_names{ii};
    eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% Enter model equations here


ai = 1;
tau=1;
d=1;
ex=1;
a=1;
jj=jjss;
g=ggss;

pi = 1;
pid = 1;
pim = 1;
pibard = 1;
pibarm = 1;

p = 1;
e = 1 / 1.2;
xd = 1.2;
xm = 1.2;

r = 1 / Betas;
pistar = 1.0001;
rstar = r*pistar;

m = Me;


eps_bc =    (1-Betab*Epsilonc) / (1-Epsilonc);
eps_bh = jj*(1-Betab*Epsilonh) / (1-Epsilonh);

eps_sc =    (1-Betas*Epsilonc) / (1-Epsilonc);
eps_sh = jj*(1-Betas*Epsilonh) / (1-Epsilonh);

eps_ec =    (1-Betae*Epsilonc) / (1-Epsilonc);


alphae = 1 - (1-Me*rstar/r)/2/(1-Mestar);

[ya, he, hs, hb, ns, nb, lambe, lambs, lambb, q, k] = solve12var(Nu, Mu, Alpha, xd, r, rstar, Me, Mestar, Mb, Betas, Betab, Betae, eps_ec, eps_sc, eps_sh, eps_bc, eps_bh, g, Omegax, Omegam, Eta, alphae, Delta);


i = Delta*k;

y = (1-Omegax)/(1-Omegam)*ya;
yd = (1-Omegam)*y;
yx = Omegax*ya;
ym = Omegam*y;

be = Me*alphae*q*he/r;
bestar = (1-alphae)*q*he*(1-(1-Mestar)*(1-alphae))/e/rstar;

% bb = Mb*q*hb/r;
bb = Me*q*hb/r;
bs = bb + be;

bsstar = (yx - e*ym)/(1-r)/e + bestar;



ws = Alpha * (1-Nu-Mu) * ya / ns / xd;
wb = (1-Alpha) * (1-Nu-Mu) * ya  / nb / xd;


ce = eps_ec/lambe;
cs = eps_sc/lambs;
cb = eps_bc/lambb;

c = ce+cs+cb;
tb = yx - e*ym;

ucb = lambb;
uce = lambe;
ucs = lambs;


uhs = eps_sh / hs;
uhb = eps_bh / hb;

mue = lambe/r*(1-Betae*r);
muestar = lambe/rstar*(1-Betae*r);
mub = lambb/r*(1-Betab*r);



auxn_d = 1/(1-Betas*Thetad)*(ucs * ya / xd);
auxd_d = 1/(1-Betas*Thetad)*(ucs * ya);

auxn_m = 1/(1-Betas*Thetam)*(ucs * ym / xm);
auxd_m = 1/(1-Betas*Thetam)*(ucs * ym);



Ws = 1/(1-Betas) * (log(1-Epsilonc) + log(cs) + jj*(log(1-Epsilonh) + log(hs)) - ns^(1+Eta)/(1+Eta));
Wb = 1/(1-Betab) * (log(1-Epsilonc) + log(cb) + jj*(log(1-Epsilonh) + log(hb)) - nb^(1+Eta)/(1+Eta));
We = 1/(1-Betae) * (log(1-Epsilonc) + log(ce));


%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = M_.endo_names{ii};
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end
