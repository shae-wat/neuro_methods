% input layer feeding k-winners-take-all, learning implemented on output 
% Foundation from HW 6's neural integration in the mutual inhibition circuit

close all; clear all;
load Dict_K_50_lambda_0.05.mat;

N=50; % num neurons
tau = 1; % time is in multiples of tau
dt = 0.02; % fraction of tau
Trun = 40; % multiples of tau
Nrun = Trun/dt; % num timesteps
time=0:dt:Trun-dt; % multiples of tau
s=zeros(N,Nrun); %full trace of synaptic activity of cells
basis_length = 1440;


%===========weights for basis layer's k winners===========

beta = 2;  
I = beta*ones(N,N);
 
alpha = 40; % steepness of sigmod nonlinearity
A = alpha*(eye(N));

W = A - I;


%===========input layer===========

num_movements = 1; %number of motor units to simulate
input =  (1);   %range = 1:num_movents, represents set of motor units
B = eye(num_movements); %allows selection of motor unit

num_rounds = 20; %number of learning rounds to perform for each movement

% TODO*graph* each col will hold the O during learning of each motor unit 
output = ones(basis_length, num_rounds);

%Summation of weighted basis functions from original data set
alpha_weights = [20, 30, 40];
summation = alpha_weights(1)*D(:,1) +  alpha_weights(2)*D(:,2) +  alpha_weights(3)*D(:,3);
desired_movement = summation + randn(basis_length, 1); %adds Gaussian noise
%d_movement_matrix = ones(basis_length, num_rounds);
learning_rounds = (1:num_rounds);

%===========Learn each motor unit===========
for mu=1:num_movements
    %selects desired movement/motor unit to learn
    this_movement = input(mu);
    activated_motor_unit = B(:, this_movement);
    %K is the  weights connecting input neurons to the mid layer
    K = rand(N,num_movements); %cols of weights learned for each motor unit

    hold on;
    plot(learning_rounds, sum(desired_movement), '--ob');
    
    
    %===========Learning rounds for each motor unit===========
    for round=1:num_rounds
    
        b = K * activated_motor_unit; %get appropriate col from K for this movement

        %===========basis and output layer network==========
        % k winners are the neurons represented in the last column of s
        % D is the given dictionary of basis functions repping sparse movement
        % O is the output influenced by K weights

        %===========minus phase===========
        for i=1:Nrun-1,

            %output is the dict of basis functions * this step of k winners
            O = D*s(:,i);
            %weighted input to neuron + influence by input layer 
            G = D'*O + W*s(:,i) + b;
            frate = zeros(N,1) + G.*(G>0 & G<1) + ones(N,1).*(G>1);

            %calculate next time step based on frate of this time step
            s(:,i+1) = s(:,i) + dt*(-s(:,i)+frate);

        end
        minus_phase_activation = s(:,end); %k winners
        %add calculated output vector to a matrix for future comparisons
        O = D*s(:,end); %learned output
        %d_movement_matrix(:,round) = desired_movement;


        %===========plus phase===========
        for i=1:Nrun-1,

            %weighted input to neuron + influence by input layer 
            G = D'*desired_movement + W*s(:,i) + b;
            frate = zeros(N,1) + G.*(G>0 & G<1) + ones(N,1).*(G>1);

            %calculate next time step based on frate of this time step
            s(:,i+1) = s(:,i) + dt*(-s(:,i)+frate);

        end
        plus_phase_activation = s(:,end); %k winners

        %should should differences
        %plot([minus_phase_activation, plus_phase_activation]);


        %===========update basis function weights===========
        dynamic_thres_param = .5;
        rev_dir_pt = .1;
        learning_rate = 0.1;
        delta = f_xcal(plus_phase_activation, minus_phase_activation, activated_motor_unit, dynamic_thres_param, rev_dir_pt);

        %update weights from input layer to basis layer
        K = K - learning_rate.*delta; 
        
        
        %difference between desired_movement and each output
        
        plot(learning_rounds, sum(desired_movement - O), '--or');
        %hold on;
        

    
    end
    

end