%% Initialization
clear; close all;

% Parameters
N = 1000; % Total number of neurons
Np = 200; % Number of neurons in the integration subnetwork
tau_m = 20; % Membrane time constant (ms)
theta = 0.1; % Spiking threshold
tau_s = 5; % Synaptic conductance time constant (ms)
T = 1000; % Total simulation time (ms)
dt = 0.1; % Time step (ms)
timesteps = T/dt; % Number of time steps

% Initialize membrane potentials and synaptic currents
x = zeros(N, timesteps); % Membrane potential for each neuron over time
p = zeros(N, timesteps); % Synaptic input for each neuron over time
r = zeros(N, timesteps); % Spiking rate for each neuron over time

% Connectivity matrices
W = rand(N, N)./sqrt(N); % Synaptic weight matrix for the whole network
W(1:Np, 1:Np) = rand(Np, Np)./sqrt(Np); % Subnetwork weight matrix

%% External input
input_duration = 100; % Duration of input pulses (ms)
pulse_isi = 20; % Inter-stimulus interval (ms)
input_times = [100, 100+input_duration+pulse_isi, ...
               100+2*(input_duration+pulse_isi), ...
               100+3*(input_duration+pulse_isi)];
input_neurons = randperm(N, round(0.25*N)); % Random 25% of neurons receive input

% Generate the smoothened step function for the external input
s = zeros(N, timesteps);
for i = 1:length(input_times)
    pulse_start = round(input_times(i)/dt);
    pulse_end = round((input_times(i) + input_duration)/dt);
    s(input_neurons, pulse_start:pulse_end) = 1;
end

% Smooth the input pulses to create a more gradual onset and offset
sigma = 100; % Smoothing factor (adjustable for different smoothness)
for i = 1:length(input_neurons)
    s(input_neurons(i), :) = smooth(s(input_neurons(i), :), sigma);
end

figure;
plot(dt:dt:T, sum(s,1));

%% Calculate the network time constant
lambda_max = max(real(eig(W))); % Calculate the largest eigenvalue of the synaptic weight matrix W
tau_n = tau_s / abs(1 - lambda_max);

%% Simulation loop
% for t = 1:timesteps
%     % Update synaptic inputs
%     p = p + dt * (-p/tau_s + r);
%     
%     % Update membrane potentials
%     x = x + dt/tau_m * (-x + W*p + s(:, t));
%     
%     % Check for spiking neurons
%     spiking_neurons = find(x > theta);
%     r(spiking_neurons) = 1; % Set spiking rate to 1
%     x(spiking_neurons) = 0; % Reset membrane potential after spike
%     
%     % Reset spiking rate to 0 for next timestep
%     r = r * 0;
% end

% Simulation loop
dtau_s = tau_s/dt;
dtau_m = tau_m/dt;
for t = 2:timesteps % Start at 2 to access the previous timestep
    % Update synaptic inputs
    p(:, t) = p(:, t-1) + dt * (-p(:, t-1)/dtau_s + r(:, t-1));
    
    % Update membrane potentials
    x(:, t) = x(:, t-1) + dt/dtau_m * (-x(:, t-1) + W*p(:, t-1) + s(:, t));
    
    % Check for spiking neurons
    spiking_neurons = find(x(:, t) > theta);
    r(spiking_neurons, t) = 1; % Set spiking rate to 1
    x(spiking_neurons, t) = 0; % Reset membrane potential after spike
end

%% Plot results
figure;
plot(dt:dt:T, sum(r, 1));
xlabel('Time (ms)');
ylabel('Total Spiking Rate');
title('Total Spiking Activity Over Time');

% Plot results: Example Neuron Membrane Potential Over Time
example_neuron = randi(N); % Change this index to view different neurons
figure;
plot(dt:dt:T, x(example_neuron, :));
xlabel('Time (ms)');
ylabel('Membrane Potential');
title(['Membrane Potential of Neuron ', num2str(example_neuron)]);