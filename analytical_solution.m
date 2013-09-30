%Using the synaptic conductance (synaptic activation) equation 
%subject to a train of action potentials occuring at given spike times
%Plot of the analytical solution for gsyn(t) over the interval [0,1] seconds

t = 1000; %1000 ms = 1 second
dt = 0.1;
nsteps = t/dt; %numbr of time steps
tau_syn = 50; %ms
const_g_syn = 0.1; %height of bump

n = 41; %number of spikes
spike_times = [2, 27, 46, 75, 87, 111, 128, 148, 170, 200, 223, 245, 263, 283, 298, 317, 344, 362, 382, 404, 425, 444, 464, 478, 504, 522, 538, 555, 574, 595, 612, 629, 649, 673, 695, 717, 740, 762, 784, 799, 818];

g_syn_init = 0;
g_syns = g_syn_init; %vector to hold evaluated g_syn(t) to plot
time = 0:dt:(t-dt); %vector of time steps

%for loop over spike times
for i=1:n
    %vector of 0s before time and 1s after time of spike
    g_set = spike_times(i) <= time;
    %bump and decay of a single spike over all time
    g = const_g_syn*exp((spike_times(i)-time)/tau_syn);
    %use g_set to set pre-spike values to 0
    %then add to summation of all spikes over time
    g_syns = g_syns + g.*g_set;
end

%(1a) Analytical solution for g_syn over interval
plot(time, g_syns);
xlabel('time');
ylabel('gsyn(t)');

