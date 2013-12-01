%Mutual integration in the mutual inhibition circuit
% linear transfer function in the circuit
% frate = W*s(:,i) + b
% b = sinusoidal function over time

close all; clear all;
N=2; % num neurons
tau = 1; % time is in multiples of tau
dt = 0.02; % fraction of tau
Trun = 100; % multiples of tau
Nrun = Trun/dt; % num timesteps
time=0:dt:Trun-dt; % multiples of tau
alpha = 5; % steepness of sigmod nonlinearity
w=1; 
W=[[0,-w];[-w,0]];
s=zeros(2,Nrun); % synaptic activity of cells full trace

% initial conditions
s(1,1) = 0.2; 
s(2,1) = 1.2;

b = sin((1:Nrun)*dt*2*pi/20)*4;

%linear transfer function
f = @(x) x;

for i=1:Nrun-1,
    
    frate= f(W*s(:,i) + b(i));
    s(:,i+1) = s(:,i) + dt*(-s(:,i)+frate);

end

figure(1);
plot(time,s');
xlabel('time');
title('Mutual inhibition circuit with sinsusoidal b_1 = b_2');

%neurons maintain their stability due to matching changes in b input
