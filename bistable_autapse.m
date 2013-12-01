%Graphical stability analysis over different values of alpha
%show how the 2 stable points grow away from the point x = 1/2

%over x
dx = .05;
x = -.5:dx:1.5;
%x curve
funct_x = x;
%f(x)
alpha = [2, 3.6, 4, 4.4, 7];
f = zeros(length(x), length(alpha));

for i=1:length(alpha)
     f(:,i) = exp(alpha(i)*(x-1/2))./(1+exp(alpha(i)*(x-1/2)));
end
figure(1);
clf
hold on
%alpha = 2 in magenta
plot(x, f(:,1), 'm');
%alpha = 3.6 in red
plot(x, f(:,2), 'r');
%alpha = 4 in green
plot(x, f(:,3), 'g');
%alpha = 4.4 in magenta
plot(x, f(:,4), 'm');
%alpha = 7 in cyan
plot(x, f(:,5), 'c');
%x
plot(x, funct_x, 'k');
xlabel('x');
title('f(x) stretching over x');
hold off
 
