Inv = zeros(64,16);

for i = 1:16
    set_param_value('Rpi',CETable(5,i));
    set_param_value('Rs',CETable(6,i));
    set_param_value('Ry',CETable(7,i));
    set_param_value('Rq',CETable(8,i));
    set_param_value('Phiq',CETable(9,i));
    set_param_value('Phiy',CETable(10,i));
    set_param_value('Phib',CETable(11,i));
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, M_.endo_names);
    Inv(:,i) = oo_.mean;
end

