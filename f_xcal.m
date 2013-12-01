function [dw] = f_xcal(plus_phase_activation, minus_phase_activation, activated_motor_unit, rev_di)



theta_d = rev_di;

%theta_p = dynamic_thres_param;



a = plus_phase_activation*activated_motor_unit';%plusphase

b = minus_phase_activation*activated_motor_unit';%minus phas





dw = zeros(size(a));



ps = a>(b.*theta_d);

dw = a-b;

dw(ps) = a(ps)-b(ps);

 dw(~ps) = -b(~ps)*(1-theta_d)/theta_d;

