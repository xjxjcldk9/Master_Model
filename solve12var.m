function [ya, he, hs, hb, ns, nb, lambe, lambs, lambb, q, k]= solve12var(Nu, Mu, Alpha, xd, r, rstar, Me, Mestar, Mb, Betas, Betab, Betae, eps_ec, eps_sc, eps_sh, eps_bc, eps_bh, g, Omegax, Omegam, Eta, alphae, Delta)
    %%%%%%%%%%%1: ya, 2: he, 3: hs, 4: hb, 5: ns, 6: nb, 7: lambe, 8: lambs, 9: lambb, 10: q
    fun = @(x)eqns(x, Nu, Mu, Alpha, xd, r, rstar, Me, Mestar, Mb, Betas, Betab, Betae, eps_ec, eps_sc, eps_sh, eps_bc, eps_bh, g, Omegax, Omegam, Eta, alphae, Delta);
    init = ones(1,11);
    init(1) = 1.5877;
    init(2) = 0.0655;
    init(3) = 0.8807;
    init(4) = 0.0537;
    init(5) = 0.9159;
    init(6) = 1.0204;
    init(7) = 4.1244;
    init(8) = 1.0342;
    init(9) = 5.8530;
    init(10) = 92.2262;
    init(11) = 7.0973;
  
    %'Display','iter' 'FunctionTolerance', 1e-10, @lsqnonlin,'PlotFcn','optimplotfval'
    %opts = optimoptions(@fsolve, 'MaxFunctionEvaluations', 300000, 'MaxIterations', 3000, 'Diagnostics', 'off');
    opts = optimoptions( @lsqnonlin, 'FunctionTolerance', 1e-11, 'Display','off');
    %X = fsolve(fun, init, opts);
    lower = zeros(1,11);
    upper =  ones(1,11)*600;
 
    X = lsqnonlin(fun, init, lower, upper, opts);
    %X = lsqnonlin(fun, init,ones(1,14)*(-200),ones(1,14)*200, opts);
    
    %disp(X);

    ya = X(1);
    he = X(2);
    hs = X(3);
    hb = X(4);
    ns = X(5);
    nb = X(6);
    lambe = X(7);
    lambs = X(8);
    lambb = X(9); 
    q = X(10);
    k = X(11);
   
end


function F = eqns(x, Nu, Mu, Alpha, xd, r, rstar, Me, Mestar, Mb, Betas, Betab, Betae, eps_ec, eps_sc, eps_sh, eps_bc, eps_bh, g, Omegax, Omegam, Eta, alphae, Delta)
    %%%%%%%%%%%1: ya, 2: he, 3: hs, 4: hb, 5: ns, 6: nb, 7: lambe, 8:
    %%%%%%%%%%%lambs, 9: lambb, 10: q, 11:k
    F(1) = x(11)^Mu *x(2)^Nu * (x(5)^(Alpha) * x(6)^(1-Alpha))^(1-Nu-Mu) - x(1);
    F(2) = (Nu+Mu)*x(1)/xd + (1-r)*(Me*alphae/r+(1-alphae)*(1-(1-Mestar)*(1-alphae))/rstar)*x(10)*x(2) - eps_ec/x(7)-Delta*x(11);
    F(3) = Me*(1-Betae*r)/r + Betae*(Nu*x(1)/x(10)/x(2)/xd + 1)-1;
    F(4) = eps_sh/x(3)/x(10) - (1-Betas)*x(8);
    F(5) = Alpha*(1-Nu-Mu)*x(1)/xd*x(8)-(x(5))^(1+Eta);
%     F(6) = (1-Alpha)*(1-Nu-Mu)*x(1)/xd + (1-r)*Mb*x(10)*x(4)/r - eps_bc/x(9);
    F(6) = (1-Alpha)*(1-Nu-Mu)*x(1)/xd + (1-r)*Me*x(10)*x(4)/r - eps_bc/x(9);
%     F(7) = eps_bh/(x(4)*x(10)) + x(9)*(1/r*(1-Betab*r)*Mb-(1-Betab));
F(7) = eps_bh/(x(4)*x(10)) + x(9)*(1/r*(1-Betab*r)*Me-(1-Betab));
    F(8) = (1-Alpha)*(1-Nu-Mu)*x(1)/xd*x(9)-(x(6))^(1+Eta);
    F(9) = eps_ec/x(7) + eps_sc/x(8) + eps_bc/x(9) + Delta*x(11) - (1-g)*(1-Omegax)/(1-Omegam)*x(1);
    F(10) = x(2) + x(3) + x(4) - 1;
    F(11) = Betae*(Mu*x(1)/xd/x(11) + (1-Delta)) - 1;
end
