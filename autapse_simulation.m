%Numerical simulation of an autapse described by a given equation
%For different values of stimulation b, steady-state x is plotted
 
tau = 1; % time is in multiples of tau
dt = 0.02; % fraction of tau
Trun = 100; % multiples of tau
Nrun = Trun/dt; % num timesteps
time=0:dt:Trun-dt; % multiples of tau
W=1; 
 
alpha = 8; % steepness of sigmod nonlinearity
 
%large to small 
b1=(10:-.1:-10); % instantaneous input to cells
%small to large 
b2 = (-10:.1:10);
  
s1=zeros(length(b1),Nrun); % synaptic activity of cells full trace
s2=zeros(length(b2),Nrun);
 
g1=zeros(1,Nrun); % instantaneous total conductance into cells
frate1=g1; % instantaneous firing rate of cells
g2=zeros(1,Nrun); % instantaneous total conductance into cells
frate2=g2; % instantaneous firing rate of cells
expfact1=0; % exp term in sigmoid transfer function
expfact2=0;
 
 
% initial conditions
s1(1,1) = 1.2;
s2(1,1) = 1.2;
 
%loop over all b's
for j=1:length(b1)
 
    if j~=1
        s1(j,1) = s1(j-1,end);
        s2(j,1) = s2(j-1,end);
    end
 
    g1(1) = W*s1(j,1)+b1(j); 
    expfact1=exp(alpha*(g1(1)-1/2)); 
    frate1 = expfact1./(1+expfact1);
    g2(1) = W*s2(j,1)+b2(j); 
    expfact2=exp(alpha*(g2(1)-1/2)); 
    frate2 = expfact2./(1+expfact2);
 
    for i=1:Nrun-1,
        s1(j,i+1) = s1(j,i) + dt*(-s1(j,i)+frate1+b1(j));
        s2(j,i+1) = s2(j,i) + dt*(-s2(j,i)+frate2+b2(j));
 
        g1(i+1) = W*s1(j,i+1)+b1(j); 
        expfact1 = exp(alpha*(g1(i+1)-1/2));
        frate1 = expfact1./(1+expfact1);
        g2(i+1) = W*s2(j,i+1)+b2(j); 
        expfact2 = exp(alpha*(g2(i+1)-1/2));
        frate2 = expfact2./(1+expfact2);
    end
 
end
 
figure(2);
clf
hold on
plot(b1,s1(:,99),'b');
plot(b2,s2(:,99),'m');
xlabel('b');
ylabel('x bar');
title('autapse numerically simulated sweeps in b');
hold off
