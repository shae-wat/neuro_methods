%Work with given dictionary D

load Dict_K_50_lambda_0.05.mat;

figure(1);
plot(D(:,50));
title('basis function representing 24 muscles over 60 steps');

%the point of our project is to make the basis functions of fewer muscles
%meaningful to movement, specifically as defining motor units
figure(2);
plot(D(1:60,50));
title('meaningless = basis function representing a muscle over 60 steps');
