% input layer feeding k-winners-take-all, learning implemented on output 
% Foundation from HW 6's neural integration in the mutual inhibition

% circuit
load Dict_K_50_lambda_0.05.mat;

%mean subtract, and normalize columns D
normOfBasisFunctions = 8;
D = bsxfun(@(a,b)a-b,D,mean(D,1));
D = bsxfun(@(a,b)a./b,D,sum(D.^2,1)) * normOfBasisFunctions; 

figure(1)
clf
hold on

N=50; % num neurons
tau = 1; % time is in multiples of tau
dt = 0.02; % fraction of tau
Trun = 40; % multiples of tau
Nrun = Trun/dt; % num timesteps
time=0:dt:Trun-dt; % multiples of tau
s=zeros(N,Nrun); %full trace of synaptic activity of cells
basis_length = 1440;


%===========weights for basis layer's k winners===========
backprojectionScale = 1/1;
learning_rate = 0.005;
beta = 0.012;  
I = beta*ones(N,N);
 
alpha = 0.6; % steepness of sigmod nonlinearity
A = alpha*(eye(N));

W = A - I;

%===========input layer===========

num_movements = 1; %number of motor units to simulate
input =  (1);   %range = 1:num_movents, represents set of motor units
B = eye(num_movements); %allows selection of motor unit

num_rounds = 45; %number of learning rounds to perform for each movement

% TODO*graph* each col will hold the O during learning of each motor unit 
output = ones(basis_length, num_rounds);

%Summation of weighted basis functions from original data set
%%alpha_weights = [20, 30, 40];
%%summation = alpha_weights(1)*D(:,1) +  alpha_weights(2)*D(:,2) +  alpha_weights(3)*D(:,3);
alphas_to_use = randperm(N);%[20, 30, 40];
alphas_to_use = alphas_to_use(1:5);
alpha_weights = [1; 1; 1; 0.5; 0.5];
desired_activity = zeros(N,1);
desired_activity(alphas_to_use) = alpha_weights;
summation = D(:,alphas_to_use)*alpha_weights;
desired_movement = summation + 0.01*randn(basis_length, 1); %adds Gaussian noise
%d_movement_matrix = ones(basis_length, num_rounds);
learning_rounds = (1:num_rounds);

timeToKeepInputOn = Nrun / 2;
%===========Learn each motor unit===========
for mu=1:num_movements
    %selects desired movement/motor unit to learn
    this_movement = input(mu);
    activated_motor_unit = B(:, this_movement);
    %K is the  weights connecting input neurons to the mid layer
    K = 0.1*rand(N,num_movements); %cols of weights learned for each motor unit
    
    %===========Learning rounds for each motor unit===========
    for round=1:num_rounds
    
        b = K * activated_motor_unit; %get appropriate col from K for this movement

        %===========basis and output layer network==========
        % k winners are the neurons represented in the last column of s
        % D is the given dictionary of basis functions repping sparse movement
        % O is the output influenced by K weights

        %===========minus phase===========
        for i=1:Nrun-1,
            if(i > timeToKeepInputOn)
                b(:) = 0;
            end
            %output is the dict of basis functions * this step of k winners
            O = D*s(:,i);
            %weighted input to neuron + influence by input layer 
            G = backprojectionScale*D'*O + W*s(:,i) + b;
            frate = 1e-3 + G.*(G>0 & G<1) + ones(N,1).*(G>1);
            
            %calculate next time step based on frate of this time step
            s(:,i+1) = s(:,i) + dt*(-1/tau*s(:,i)+frate);

        end
        minus_phase_activation = s(:,end); %k winners
        %add calculated output vector to a matrix for future comparisons
        O = D*s(:,end); %learned output
        %d_movement_matrix(:,round) = desired_movement;


        %===========plus phase===========
        for i=1:Nrun-1,
            if(i > timeToKeepInputOn)
                b(:) = 0;
            end
            %weighted input to neuron + influence by input layer 
            G = backprojectionScale*D'*desired_movement + W*s(:,i) + b;
            frate = 1e-3 + G.*(G>0 & G<1) + ones(N,1).*(G>1);

            %calculate next time step based on frate of this time step
            s(:,i+1) = s(:,i) + dt*(-1/tau*s(:,i)+frate);

        end
        plus_phase_activation = s(:,end); %k winners

        %should should differences
        %plot([minus_phase_activation, plus_phase_activation]);


        %===========update basis function weights===========
        dynamic_thres_param = .5;
        rev_dir_pt = 0.2;
        delta = f_xcal(plus_phase_activation, minus_phase_activation, activated_motor_unit, rev_dir_pt);

        %update weights from input layer to basis layer
        K = K + learning_rate.*delta; 
        
        
        %difference between desired_movement and each output
        figure(1)
        plot(round, sum((desired_movement - O).^2), 'x');
        title('Difference between learned movement and desired movement');
        xlabel('rounds');
        ylabel('distance from desired movement');
        figure(2)
        plot([desired_activity,  minus_phase_activation, plus_phase_activation]);
        title('Plus and minus learning phase activations and target activity');
        legend({'Target activity','Minus activation', 'Plus activation'},'Location','NorthWest');
        xlabel('basis functions');
        ylabel('alpha weight values');
        figure(3)
        plot(plus_phase_activation-minus_phase_activation);
        title('Difference between plus-phase and minus phase basis fxns');
        xlabel('basis functions');
        ylabel('differences in alpha activation values');
        pause(0.5);
    end
end
