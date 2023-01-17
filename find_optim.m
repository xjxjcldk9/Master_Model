%%Give a table for all optimal parameters and welfare

load bigtable.mat

for i = 5:6
    assign = int2bit(i,4);
    [wws, wwb, wwe, x, val] = calwel(assign(1), assign(2), assign(3), assign(4));
    BigTable(:,i+1) = [wws wwb wwe x(1:end) val] ;
    disp(i);
end

save bigtable.mat




function [wws, wwb, wwe, x, val] = calwel(Phiq, Phiy, Phib, Rq)
    global oo_ M_ options_
    
    
    set_param_value('Me',0.85);
    
 
      
    
    options_.noprint=1; %關掉通知
    options_.order=2;
    options_.irf=0;
    
    
    ub = 20;
    lb = 0;

    
        

    Phiq_lb=0;
    Phiq_ub=0;
    
    Phiy_lb=0;
    Phiy_ub=0;
    
    Phib_lb=0;
    Phib_ub=0;
    
    Rq_lb=0;
    Rq_ub=0;
    

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
    if Rq == 1
        Rq_ub=ub;
        Rq_lb=lb;
    end

    
    [x,val, flgpso, opso] = particleswarm(@loss, 7, [lb,-ub,lb,Phiq_lb,Phiy_lb,Phib_lb,Rq_lb],[ub,ub,ub,Phiq_ub,Phiy_ub,Phib_ub,Rq_ub]);
    
    
    %要evaluate多少次
    disp(opso.funccount);

    set_param_value('Rpi',x(1));
    set_param_value('Rs', x(2));
    set_param_value('Ry', x(3));
    
    
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

    set_param_value('Phiq',xopt(4));
    set_param_value('Phiy',xopt(5));
    set_param_value('Phib',xopt(6));
    set_param_value('Rq',xopt(7));

    
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, {'pi', 'e', 'y', 'bs', 'q'});
    L = sum(sqrt(diag(oo_.var)) ./ abs(oo_.mean));


end



