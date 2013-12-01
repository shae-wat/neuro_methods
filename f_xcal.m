function [ f_xcal ] = f_xcal( plus_phase_activation, minus_phase_activation, activated_motor_unit, dynamic_thres_param, rev_dir_pt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% delta_{i,j} = f_xcal(plus_phase_activation(i)*activated_motor_unit(j), minus_phase_activation(i)*activated_motor_unit(j))
xy_plus = plus_phase_activation * activated_motor_unit';
xy_minus =  minus_phase_activation * activated_motor_unit';
f_xcal = zeros(50,1);
M = xy_plus.*xy_minus;

%TODO: matrix shit
f_xcal(M > dynamic_thres_param*rev_dir_pt) = M(M > dynamic_thres_param*rev_dir_pt) - dynamic_thres_param;
f_xcal(M <= dynamic_thres_param*rev_dir_pt) = -M(M <= dynamic_thres_param*rev_dir_pt)*(1-rev_dir_pt)/rev_dir_pt;


% for m=1:50
%     if (M > dynamic_thres_param*rev_dir_pt)
%         f_xcal(m) = M - dynamic_thres_param
%     else
%         f_xcal(m) = -M*(1-rev_dir_pt)/rev_dir_pt
%     end
% end






