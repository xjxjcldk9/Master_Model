%%Give a table for all optimal parameters and welfare

%load bigtable.mat

%央行的參數可以調大，每個都要乘以0.0849，最滿20也就1.7多而已
%那MP最多也就2好了

BigTable = zeros(11,16);

for i = 0:15
    assign = int2bit(i,4);
    [wws, wwb, wwe, x, val] = calwel(assign(1), assign(2), assign(3), assign(4));
    BigTable(:,i+1) = [wws wwb wwe x(1:end) val] ;
    disp(i);
end


%save bigtable.mat




function [wws, wwb, wwe, x, val] = calwel(Rq, Phiq, Phiy, Phib)
    global oo_ M_ options_
    
    
    set_param_value('Me',0.90);
    
 
    options_.noprint=1; %關掉通知
    options_.order=2;
    options_.irf=0;
    
    
    cub=20;


    ub = 2;
    lb = 0;
   
    Rq_lb=0;
    Rq_ub=0;

    Phiq_lb=0;
    Phiq_ub=0;
    
    Phiy_lb=0;
    Phiy_ub=0;
    
    Phib_lb=0;
    Phib_ub=0;
    
    
    if Rq == 1
        Rq_ub=cub;
        Rq_lb=lb;
    end

    if Phiq == 1
        Phiq_ub=ub;
        Phiq_lb=lb;
    end
    if Phiy == 1
        Phiy_ub=ub;
        Phiy_lb=lb;
    end
    if Phib == 1
        Phib_ub=ub;
        Phib_lb=lb;
    end
    

    
    [x,val, flgpso, opso] = particleswarm(@loss, 7, [lb,-cub,lb,Rq_lb, Phiq_lb, Phiy_lb, Phib_lb], [cub, cub, cub, Rq_ub, Phiq_ub, Phiy_ub, Phib_ub]);
    
    
    %要evaluate多少次
    disp(opso.funccount);

    set_param_value('Rpi',x(1));
    set_param_value('Rs', x(2));
    set_param_value('Ry', x(3));
    
    set_param_value('Rq', x(4));
    set_param_value('Phiq', x(5));
    set_param_value('Phiy', x(6));
    set_param_value('Phib', x(7));
    
    
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, {'Ws', 'Wb', 'We'});
     
    wws = oo_.mean(1);
    wwb = oo_.mean(2);
    wwe = oo_.mean(3);

end


function L = loss(xopt)
    global oo_ M_ options_

    set_param_value('Rpi',xopt(1));
    set_param_value('Rs',xopt(2));
    set_param_value('Ry',xopt(3));
    set_param_value('Rq',xopt(4));
    set_param_value('Phiq',xopt(5));
    set_param_value('Phiy',xopt(6));
    set_param_value('Phib',xopt(7));
    
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, {'pi', 'e', 'y', 'bs', 'q'});
    L = sum(sqrt(diag(oo_.var)) ./ abs(oo_.mean));


end



