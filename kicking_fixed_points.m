close all; clear all;
 
tau = 1; % time is in multiples of tau
dt = 0.02; % fraction of tau
Trun = 100; % multiples of tau
Nrun = Trun/dt; % num timesteps
time=0:dt:Trun-dt; % multiples of tau
 
alpha = 5; % steepness of sigmod nonlinearity
 
W=1; 
s=zeros(1,Nrun); % synaptic activity of cells full trace
g=zeros(1,Nrun); % instantaneous total conductance into cells
frate=g; % instantaneous firing rate of cells
 
b=0; % instantaneous input to cells
 
expfact=0; % exp term in sigmoid transfer function
 
% initial conditions
s(1) = 1/2; %initialized for the unstable point
g(1) = W*s(1)+b; 
expfact=exp(alpha*(g(1)-1/2)); 
frate = expfact./(1+expfact);     
 
% some diagnostics
figure(1); 
u=-2:.05:2; 
plot(u,exp(alpha*(u-1/2))./(1+exp(alpha*(u-1/2)))); 
 
for i=1:Nrun-1,
 
    s(i+1) = s(i) + dt*(-s(i)+frate);
 
    % halfway through, make a switch: kick the switches in opposite directions. 
    % kick duration: 5 tau
    if (i<(Trun/2/dt+5/dt))&(i>=(Trun/2/dt)),
        deltab = 0.05; 
        b = deltab;
        %b = -deltab; 
    else
        b=0; 
    end
  
    g(i+1) = W*s(i+1)+b; 
    expfact = exp(alpha*(g(i+1)-1/2));
    frate = expfact./(1+expfact);     
    
end
 
figure(3);
plot(time,s); 
xlabel('time');
ylabel('s');
title('unstable fixed point with positive kick');
