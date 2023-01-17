//****************************************************************************
//Define variables
//****************************************************************************

var 
    y           (long_name='Final output')
    ya          (long_name='homogeneous produced by entrepreneur')
    yd          (long_name='Domestic output')
    yx          (long_name='Export')
    ym          (long_name='Import')
    p           (long_name='Final real price')
    e           (long_name='Real exchange')
    xd          (long_name='domestic markup')
    xm          (long_name='foreign markup')   
            
    i           (long_name = 'investment')
    k           (long_name = 'capital')

    pi          (long_name='final good inflation')
    pid         (long_name='domestic price inflation') 
    pim         (long_name='import price inflation') 
    pibard       (long_name='Optimal domestic new price')    
    pibarm       (long_name='Optimal import new price') 
    pistar      (long_name='foreign inflation')
    
    alphae      (long_name='domestic collateral for entrepreneur')


    r           (long_name='interest rate')
    rstar       (long_name='Foreign interest rate')
    
    q           (long_name='Real house price')
    ce          (long_name='entrepreneur consumption')
    cs          (long_name='saver consumption')
    cb          (long_name='borrower consumption')
    c           (long_name='consumption')
    
    tb          (long_name='trade balance')

    he          (long_name='entrepreneur housing')
    hs          (long_name='saver housing')
    hb          (long_name='borrower housing')
    ns          (long_name='saver working hour')
    ws          (long_name='saver wage rate')
    nb          (long_name='borrower working hour')
    wb          (long_name='borrower wage rate')
    be          (long_name='entrepreneur domestic bond')
    bestar      (long_name='entrepreneur foreign bond')
    bs          (long_name='saver domestic bond')
    bsstar      (long_name='saver foreign bond')
    bb          (long_name='borrower domestic bond')

    lambe       (long_name = 'e Lagrange lamb')
    mue         (long_name = 'e Lagrange mu')
    muestar     (long_name = 'e Lagrange mustar')
    lambs       (long_name = 'saver Lagrange lamb')
    lambb       (long_name = 'b Lagrange lamb')
    mub         (long_name = 'b Lagrange mu')
    
    auxn_m  (long_name = 'Calvo import aux num')
    auxd_m   (long_name = 'Calvo import aux denom')
    auxn_d  (long_name = 'Calvo domestic aux num')
    auxd_d   (long_name = 'Calvo domestic aux denom')
    
    ai
    tau         (long_name='labor shock')
    jj          (long_name='housing shock')
    d           (long_name='intertemporal shock')
    ex          (long_name='foreign demand shock')
    a           (long_name='technology shock')
    g           (long_name = 'gov shock')
   
    uce
    ucs
    uhs
    ucb
    uhb

    m           (long_name='LTV')


    Ws          (long_name = 'Welfare saver')
    Wb          (long_name = 'Welfare borrower')
    We          (long_name = 'Welfare entre')
    ;




varexo 
    eps_a 
    eps_ai
    eps_jj 
    eps_d 
    eps_tau 
    eps_g 
    eps_pistar
    eps_rstar 
    eps_ex 
    eps_u 
;
    
parameters 
    Kappaz  (long_name='domestic and foreign good elasiticty')
    Kappax  (long_name='domestic and import good elasiticty')
    Omegax  (long_name='long run export proportion')
    Omegam  (long_name='long run import proportion')
    Xi      (long_name='Domestic elasticity')
    Thetad     (long_name='Domestc Fixed price')
    Thetam     (long_name='Import Fixed price')
    Betas       (long_name='Saver discount')
    Betab       (long_name='Borrower discount')
    Betae       (long_name='e discount')
    Mu
    Nu      (long_name='Housing share')
    Alpha      (long_name='Saver work share')
    Epsilonc      (long_name='Consumption habit')
    Epsilonh      (long_name='Housing habit')
    Eta      (long_name='Labor supply elasticity')
    Delta      (long_name='Depreciation rate')
    Me      (long_name='Domestic LTV ')
    Mestar      (long_name='Foreign LTV ')
    Mb
    Phi     (long_name = 'capital adjustment cost')
    Phistar     (long_name = 'foreign bond adjustment cost')
    RR          (long_name = 'Taylor R')
    Rpi          (long_name = 'Taylor pi')
    Rs          (long_name = 'Taylor s')
    Ry          (long_name = 'Taylor y')
    Phib
    Rq
    Phiq
    Phiy


    Rhoa
    Rhod
    Rhoai
    Rhojj
    jjss
    ggss
  
    Rhotau
    Rhog
    Rhorstar
    Rhoex
    Rhopistar
    ;
//****************************************************************************
//Set parameter values
//****************************************************************************

Kappaz  =   1.6220;
Kappax  =   3.4054;
Omegax  =   0.6014;
Omegam  =   0.6582;
Xi      =   6;
Thetad  =   0.4742;
Thetam  =   0.5367;
Betas   =   0.9978;
Betab   =   0.97;
Betae   =   0.97;
Mu      =   0.3;
Nu      =   0.03;
Alpha   =   0.7926;
%Alpha   =   0.01;
%Epsilonc =  0.0;
Epsilonc =  0.4216;
Epsilonh =  0.0;
Eta      =  2.6352;
Delta    =  0.025;
Me       =  0.85;
Mestar   =  0.5;
Mb       =  0.85;
Phi      =  13.4387;
Phistar  =   4.1898 * 1e-5;
RR       =  0.9151;
Rpi      =  1.1955;
Rs       =  -0.0968;
Ry       =  0.1707;
jjss    =   0.1848;
ggss    =   0.1530;


Rq = 0;
Phiq = 0;
Phiy = 0;
Phib = 0;

Rhoa        = 0.8301;    
Rhod        = 0.8712;
Rhoai       = 0.9019;
Rhojj       = 0.9933;
Rhotau      = 0.8809;
Rhog        = 0.9233;
Rhorstar    = 0.8085;
Rhoex       = 0.8551;
Rhopistar   = 0.4099;






//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************

model;
    [name = 'Welfare saver recursive']
    Ws = d*(log(cs-Epsilonc*cs(-1)) + jj*log(hs-Epsilonh*hs(-1)) - tau/(1+Eta)*ns^(1+Eta)) + Betas*Ws(+1);
    [name = 'Welfare borrower recursive']
    Wb = d*(log(cb-Epsilonc*cb(-1)) + jj*log(hb-Epsilonh*hb(-1)) - tau/(1+Eta)*nb^(1+Eta)) + Betab*Wb(+1);
    [name = 'Welfare entre recursive']
    We = d*(log(ce-Epsilonc*ce(-1))) + Betae*We(+1);


    uce = d / (ce - Epsilonc*ce(-1)) - Betae*d(+1)*Epsilonc / (ce(+1) - Epsilonc*ce);
    ucs = d / (cs - Epsilonc*cs(-1)) - Betas*d(+1)*Epsilonc / (cs(+1) - Epsilonc*cs);
    uhs = jj*d / (hs - Epsilonh*hs(-1)) - Betas*jj(+1)*d(+1)*Epsilonh / (hs(+1) - Epsilonh*hs);
    ucb = d / (cb - Epsilonc*cb(-1)) - Betab*d(+1)*Epsilonc / (cb(+1) - Epsilonc*cb);
    uhb = jj*d / (hb - Epsilonh*hb(-1)) - Betab*jj(+1)*d(+1)*Epsilonh / (hb(+1) - Epsilonh*hb);



    [name = 'Domestic good demand']
    yd = (1-Omegam) * p^(Kappaz) * y;
    [name = 'import good demand']
    ym = (Omegam) * (p/(xm*e))^(Kappaz) * y;
    
    [name = 'index']
    p^(1-Kappaz) = (1-Omegam) + (Omegam) * (xm*e)^(1-Kappaz);

    [name = 'domestic inflation']
    pi/pid = p/p(-1);

    [name = 'import inflation']
    pim/pid = xm*e/xm(-1)/e(-1);

    [name = 'Export Demand']
    yx = Omegax * STEADY_STATE(e^(-Kappax)*ya) * (e)^(Kappax) * ex ;
    
    [name = 'wholesale equilibrium']
    ya = yd + yx;
    

    [name = 'Calvo wholesale']
    Xi * auxn_d = (Xi-1) * auxd_d;
    [name = 'auxn_d']
    auxn_d = ucs * ya / xd + Betas*Thetad*pid(+1)^(Xi)*auxn_d(+1);
    [name = 'auxd_d']
    auxd_d = ucs * pibard * ya + Betas*Thetad*(pid(+1))^(Xi-1)*pibard/pibard(+1)*auxd_d(+1);
    
    [name = 'domestic price evolution']
    1 = Thetad * pid^(Xi-1) + (1-Thetad)*(pibard)^(1-Xi);

    
    [name = 'Calvo import']
     Xi * auxn_m = (Xi-1) *  auxd_m;
    [name = 'auxn_m']
    auxn_m = ucs * ym / xm + Betas*Thetam*pim(+1)^(Xi)*auxn_m(+1);
    [name = 'auxd_m']
    auxd_m = ucs * pibarm * ym + Betas*Thetam*(pim(+1))^(Xi-1)*pibarm/pibarm(+1)*auxd_m(+1);
    
    [name = 'Import price evolution']
    1 = Thetam * pim^(Xi-1) + (1-Thetam)*(pibarm)^(1-Xi);
    
   
    [name = 'Production function']
    ya = a * (k(-1))^(Mu)*(he(-1))^(Nu) * ((ns)^(Alpha)*(nb)^(1-Alpha))^(1-Nu-Mu);

    [name = 'e budget']
    p*(ce + i) + q*(he-he(-1)) + ws*ns + wb*nb = ya/xd + (be - r(-1)/pid*be(-1)) + e*(bestar - rstar(-1)/pistar*bestar(-1) - Phistar/2*(bestar - STEADY_STATE(bestar))^2);
    
    [name = 'EOM k']
    k = ai*i + (1-Delta)*k(-1) - Phi/2*(k/k(-1)-1)^2*k(-1);

    [name = 'e domestic collat']
    r*be = m * alphae * pid(+1) * q(+1) * he;
    
    [ame = 'e foreign collat']
    e*rstar*bestar = pid(+1)*(1-alphae)*q(+1)*he*(1-(1-Mestar)*(1-alphae)*q(+1)*he/STEADY_STATE(q*he));

    [name = 'e FOC c']
    lambe * p = uce;

    [name = 'e FOC h']
    lambe*q = Betae*lambe(+1)*(Nu*ya(+1)/he/xd(+1) + q(+1)) + m*mue*alphae*q(+1)*pid(+1) + muestar*(1-alphae)*pid(+1)*q(+1)*(1-2*(1-Mestar)*(1-alphae)*q(+1)*he/STEADY_STATE(q*he));

    [name = 'e FOC doemstic bond']
    lambe = Betae * lambe(+1) * r / pid(+1) + mue*r;

    [name = 'e FOC foreign bond']
    (1-Phistar*(bestar - STEADY_STATE(bestar)))*lambe*e = Betae*lambe(+1)*e(+1)*rstar/pistar(+1) + muestar*e*rstar;

    [name = 'e FOC saver hour']
    ws = Alpha * (1-Nu-Mu) * ya / ns/xd;

    [name = 'e FOC borrower hour']
    wb = (1-Alpha) * (1-Nu-Mu) * ya / nb/xd;

    [name = 'e FOC k']
    p/ai*lambe*(1+Phi*(k-k(-1))) = Betae*lambe(+1)*(Mu*ya(+1)/xd(+1)/k + p(+1)/ai(+1) *( (1-Delta)-Phi/2*(1-(k(+1)/k)^2)));

    [name = 'e FOC alphae']
    mue*m*pid(+1)*q(+1) = muestar*pid(+1)*q(+1)*(1-2*(1-Mestar)*(1-alphae)*q(+1)*he/STEADY_STATE(q*he));
    
    [name = 's FOC c']
    lambs * p = ucs;
    [name = 's FOC h']
    lambs * q = Betas*lambs(+1)*q(+1) + uhs;
 
    [name = 's FOC hour']
    d*tau*(ns)^Eta = lambs*ws;
     
    [name = 's domestic bond']
    lambs = Betas*lambs(+1)*r/pid(+1);

    [name = 's foreign bond']
    (1 + Phistar*(bsstar - STEADY_STATE(bsstar)))*lambs*e = Betas*lambs(+1)*e(+1)*rstar/pistar(+1);

    [name = 'b bc']
    p*cb + q*(hb-hb(-1)) = wb*nb + (bb - r(-1)/pid*bb(-1));

    [name = 'b collat']
    r*bb = m*q(+1)*pid(+1)*hb;
   
    [name = 'b FOC c']
    lambb * p = ucb;

    [name = 'b FOC h']
    lambb * q = Betab*lambb(+1)*q(+1) + mub*m*q(+1)*pid(+1) + uhb;
   
    [name = 'b FOC hour']
    d*tau*(nb)^Eta = lambb*wb;

    [name = 'b domestic bond']
    lambb = Betab*lambb(+1)*r/pid(+1) + mub*r;

    [name = 'Taylor']%要加不同regime
   
    log(r) = RR*log(r(-1)) + (1-RR) * (Rpi*log(pi) + Ry *log(y/y(-1)) +  Rs*log(e/e(-1)) + Rq*log(q/q(-1)) + log(STEADY_STATE(r)))  + eps_u;
 
    [name = 'Market Clear']
    p * y = p*(cs+cb+ce + i) + e*(Phistar/2*(bestar - STEADY_STATE(bestar))^2 + Phistar/2*(bsstar - STEADY_STATE(bsstar))^2) + p*g*y;

    [name = 'housing clear']%43
    1 = hs + hb + he;

    [name = 'bond clear']%44
    bs = bb + be;

    [name = 'balanced']%45
    e * ((bsstar - rstar/pistar*bsstar(-1)) + (bestar - rstar/pistar*bestar(-1))) = yx - e*ym;
  
    log(a) = Rhoa*log(a(-1)) + eps_a;%55
    log(ai) = Rhoai*log(ai(-1)) + eps_ai;

    log(jj) = (1-Rhojj)*log(jjss) + Rhojj*log(jj(-1)) + eps_jj;%56
    log(d) = Rhod*log(d(-1)) + eps_d;%57
    log(tau) = Rhotau*log(tau(-1)) + eps_tau;%58
    
    log(g) = (1-Rhog)*log(ggss) + Rhog*log(g(-1)) + eps_g;%60    
    log(pistar) = (1-Rhopistar)*log(STEADY_STATE(pistar))+ Rhopistar*log(pistar(-1)) + eps_pistar;%61
    log(rstar) =  (1-Rhorstar)*log(STEADY_STATE(rstar)) + Rhorstar*log(rstar(-1)) + eps_rstar;%62
    log(ex) = Rhoex*log(ex(-1)) + eps_ex;%63
    

    c = cs+cb+ce;
    tb = yx-e*ym;
    


    %LTV要加盯不同東西
    m = Me*(bs/bs(-1))^(-Phib)*(q/q(-1))^(-Phiq)*(y/y(-1))^(-Phiy);

end;


shocks;
var eps_a; stderr 0.0218;
var eps_jj; stderr 0.0795; 
var eps_d; stderr 0.00321; 
var eps_tau; stderr 0.1264; 
var eps_g; stderr 0.0847; 
var eps_pistar; stderr 0.0025;
var eps_rstar ; stderr 0.0012;
var eps_ex ; stderr 0.0370;
var eps_u ; stderr 0.0008;
end;




steady;
%check ;
%model_diagnostics ;












