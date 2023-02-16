


Ws_base = BigTable(1,1);
Wb_base = BigTable(2,1);
We_base = BigTable(3,1);
Wtot_base = (1-Betas)*Ws_base + (1-Betab)*Wb_base + (1-Betae)*We_base;

CETable = zeros(12,16);

count=1;
for i = 1:16
    Ws = BigTable(1,i);
    Wb = BigTable(2,i);
    We = BigTable(3,i);
    Wtot = (1-Betas)*Ws + (1-Betab)*Wb + (1-Betae)*We;
    
    CEs = exp((1-Betas)*(Ws - Ws_base)) - 1;
    CEb = exp((1-Betab)*(Wb - Wb_base)) - 1;
    CEe = exp((1-Betae)*(We - We_base)) - 1;
    CEtot = exp((Wtot - Wtot_base)/3) - 1;
    
    CETable(1:4,count) = [CEtot CEs CEb CEe];
    CETable(5:end,count) = BigTable(4:end,i);
    
    count = count + 1;
end
