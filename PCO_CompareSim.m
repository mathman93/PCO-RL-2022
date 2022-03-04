% PCO Simulation Comparison
% Based on PCO_Sim031020
% Compares performance of PCO synchronization PRFs/algorithms
% Last Modified: 10/28/2021

%% Simulation Parameters
cycles = 201; % Maximum number of evolution cycles
alpha = 0.50; % PCO coupling strength - (0,1]
threshold = 2*pi; % Oscillator phase threshold
w_0 = 2*pi; % Nominal oscillator evolution frequency
D = 0.0; % Refractory Period
N = 6; % Number of PCOs (integer >= 2)
A_option = 'all';
A = GetTopology(N,A_option);
%A = GetTopology(N,'custom',[2:5],true); % Custom topology
eps_sync = 0.001*threshold; % epsilon-synchronization threshold

%theta_start = theta; % Store initial oscillator states (for additional simulations)

area_1 = (0.5*threshold)^2;
c = [0.07 + 0.43*(exp(-0.23*((5)-1))), 0.77 - 0.27*(exp(-0.13*((5)-1)))];
area_2 = 0.5*threshold^2 * (c(1)^2 + (1-c(2))^2 + (1-c(1))*c(1) - c(1)*(1-c(2))*(1-c(2))/(1-c(1)));

global epsilon beta; % MS parameters
epsilon = 0.01;
beta = 3;
% Estimate value of epsilon that gives same "area" as GeneralPRF
%
p_step = threshold/1000;
p = 0:p_step:threshold;
stateMS = log(1+(exp(beta)-1)*(p./threshold))./beta;
for i=1:10
    phaseMS = threshold*(exp(beta*(stateMS + epsilon))-1)./(exp(beta)-1);
    p_plus = min((phaseMS - p),(threshold-p));
    area_3 = p_step * (sum(p_plus) - 0.5*(p_plus(1)+p_plus(end)));
    epsilon = epsilon*(area_2*(1*alpha)/area_3);
end
%}

%% Simulation
rng('default');
sample_size = 100; % Number of simulations to run/compare
tau_min = 0.02; % minimum pulse delay (fraction of cycle)
tau_max = 0.04; % maximum pulse delay (fraction of cycle)
L = zeros(sample_size,1); % Initial containing arc of oscillator for each set of simulations
maxTime = zeros(sample_size,3); % Time to synchronize for each simulation
sample = 0; % Sample number index
start_time = tic; % For timing the simulations
while sample < sample_size
    omega = w_0 * (ones(1,N) + (0.0)*(rand(1,N)-0.5)); % vector of oscillator frequencies
    rand_delay = tau_min*ones(N,N) + (tau_max-tau_min)*rand(N,N); % Random matrix
    pulse_delay = 0.5 * (rand_delay + rand_delay'); % (symmetric) matrix of pulse propogation delays
    theta = rand(1,N)*threshold; % vector of (unsorted) initial node phases [0,1]
    % Run simulation and return results
    Lambda_0 = ContainingArc(theta, threshold); % Initial containing arc
    %
    if Lambda_0 < 0.5*threshold
        continue % Don't include results for "small" initial thresholds
    end
    %}
    sample = sample + 1; % Increment sample counter
    L(sample) = Lambda_0;
    %
    [time1, state1, C_Arc1] = PCO_Sync_Sim(@StandardPRF, theta, w_0, omega, cycles, A,...
                                            alpha, D, threshold, pulse_delay, eps_sync);
    [time2, state2, C_Arc2] = PCO_Sync_Sim(@GeneralPRF, theta, w_0, omega, cycles, A,...
                                            alpha, D, threshold, pulse_delay, eps_sync);
    %}
    %
    [time3, state3, C_Arc3] = PCO_Sync_Sim(@MSPRF, theta, w_0, omega, cycles, A,...
                                            alpha, D, threshold, pulse_delay, eps_sync);
    %}
    maxTime(sample,:) = [time1(end) time2(end) time3(end)];
    % Calculate percentage complete and display
    DisplayPercentageComplete(sample, sample_size, toc(start_time), 5);
end
sync_bool = maxTime < cycles;
[~,min_index] = min(maxTime,[],2);
%disp(maxTime)
disp([median(maxTime(sync_bool(:,1),1)), median(maxTime(sync_bool(:,2),2)), median(maxTime(sync_bool(:,3),3))])
disp([mean(maxTime(sync_bool(:,1),1)), mean(maxTime(sync_bool(:,2),2)), mean(maxTime(sync_bool(:,3),3))])
disp([std(maxTime(sync_bool(:,1),1)), std(maxTime(sync_bool(:,2),2)), std(maxTime(sync_bool(:,3),3))])
disp(sum(sync_bool))
disp([sum(min_index == 1), sum(min_index == 2), sum(min_index == 3)])

%% Plots
% Figure 1: PCO phases over time
%{
fig1 = figure(1);
clf
set(fig1, 'Position', [10, 300, 500, 340])
hold on
%plot(time1,state1, 'Linewidth', 1.5)
plot(time2,state2, 'Linewidth', 1.5)
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('PCO Phases (\theta)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,cycles*(threshold/w_0),0,threshold])
%}

% Figure 2: Containing Arc over time for last sample
%
fig2 = figure(2);
clf
set(fig2, 'Position', [615, 400, 500, 240])
hold on
plot(time1,C_Arc1, 'k-', 'Linewidth', 1.5) % StandardPRF
plot(time2,C_Arc2, 'r-', 'Linewidth', 1.5) % GeneralPRF
plot(time3,C_Arc3, 'b-', 'Linewidth', 1.5) % MSPRF
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('Containing Arc (\Lambda)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,max(maxTime(end,:)),0,threshold*(N-1)/N])
%set(gca, 'yscale','log')
%}

% Figure 3: Initial Containing Arc vs. time to sync
%
fig3 = figure(3);
clf
set(fig3, 'Position', [615, 50, 500, 440])
hold on
plot(L,maxTime(:,1), 'kx')
plot(L,maxTime(:,2), 'r.')
plot(L,maxTime(:,3), 'b+')
hold off
grid on
xlabelh = xlabel('Initial Containing Arc (\Lambda)');
ylabelh = ylabel('Cycles to \epsilon-synchronization');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Mirrollo-Strogatz');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northwest')
axis([min(L),threshold*(N-1)/N,0,cycles-1])
set(gca, 'yscale','log')
%}

%% Hardcoded Data Collection
alpha_values = [0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005];
alpha_values_alt = [2.0, 1.7, 1.5, 1.2, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005];
cycle_avg_MS_alt = [5.84, 7.61, 9.14, 12.45, 15.69, 17.85, 20.53, 23.95, 28.53, 35.00, 44.49, 60.56, 92.36, 97.62, 103.19, 109.39, 116.46, 124.48, 133.57, 144.22, 156.48, 170.88, 188.19, 209.55, 236.08, 270.08, 315.94, 379.91, 475.77, 635.63, 954.92, 1913.8, 4567];
cycle_avg_Opt = [5.32, 5.44, 5.54, 6.31, 6.57, 7.01, 7.75, 8.54, 9.05, 10.45, 11.47, 12.94, 14.67, 17.60, 20.99, 24.33, 31.17, 43.00, 68.11, 142.12, 303.69];
cycle_avg_RL = [6.16, 6.52, 6.72, 7.29, 7.78, 8.28, 9.09, 9.87, 11.08, 12.23, 13.93, 16.11, 17.99, 20.95, 25.14, 30.46, 37.85, 51.31, 77.11, 156.79, 311.33];
cycle_sigma_MS_alt = [3.16, 4.15, 5.0, 6.83, 8.68, 9.89, 11.42, 13.30, 15.81, 19.52, 24.69, 33.74, 51.15, 54.86, 57.77, 60.79, 65.09, 69.41, 74.34, 80.68, 87.09, 94.95, 104.48, 116.54, 130.80, 148.95, 175.26, 211.65, 265.78, 355.33, 533.24, 1071.0, 2345];
cycle_sigma_Opt = [0.51, 0.44,  0.64, 0.81, 1.19, 1.51, 2.18, 3.34, 3.56, 5.52, 5.03, 6.88, 6.48, 11.68, 13.40, 12.34, 14.46, 22.82, 41.81, 75.32, 167.58];
cycle_sigma_RL = [1.25, 1.23, 1.29, 1.34, 1.61, 2.47, 3.93, 5.25, 7.29, 8.85, 11.81, 15.27, 17.12, 20.37, 26.57, 32.28, 39.64, 55.09, 83.24, 174.95, 344.68];
sync_prob_MS_alt = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
sync_prob_Opt = [1, 0.9994, 0.9978, 0.9963, 0.9946, 0.9935, 0.9919, 0.9915, 0.9884, 0.9912, 0.9837, 0.9817, 0.9750, 0.9782, 0.9692, 0.9583, 0.9522, 0.9492, 0.9465, 0.9389, 0.9376];
sync_prob_RL = [1, 1, 1, 1, 0.9999, 0.9995, 0.9993, 0.9986, 0.9977, 0.9935, 0.9936, 0.9938, 0.9895, 0.9878, 0.9891, 0.9886, 0.9865, 0.9892, 0.9894, 0.9913, 0.9903];

% Figure 4: Probability of Sync
%{
fig4 = figure(4);
clf
set(fig4, 'Position', [115, 50, 450, 380])
%yyaxis left % Probabilities
hold on
plot(alpha_values,sync_prob_Opt, 'kx-', 'Linewidth', 1.5)
plot(alpha_values,sync_prob_RL, 'r*-', 'Linewidth', 1.5)
plot(alpha_values_alt,sync_prob_MS_alt, 'b+-', 'Linewidth', 1.5)
hold off
grid on
axis([0.0,0.2, 0.93, 1.0])
xlabelh = xlabel('Coupling Strength, l');
ylabelh_L = ylabel('Probability of Synchronization');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh_L,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Mirollo-Strogatz');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'southeast')
%axis([min(L),threshold*(N-1)/N,0,cycles-1])
%set(gca, 'yscale','log')
%set(gca, 'xscale','log')
%}

% Figure 5: Average time to epsilon-sync
%{
fig5 = figure(5);
clf
set(fig5, 'Position', [615, 50, 450, 380])
%yyaxis right % Times to Sync
hold on
plot(alpha_values, cycle_avg_Opt, 'kx-', 'Linewidth', 1.5)
plot(alpha_values, cycle_avg_RL, 'r*-', 'Linewidth', 1.5)
plot(alpha_values_alt, cycle_avg_MS_alt, 'b+-', 'Linewidth', 1.5)
hold off
grid on
axis([0.0,0.2, 4, 1000])
xlabelh = xlabel('Coupling Strength, l');
ylabelh_R = ylabel('Average Time to \epsilon-synchronization');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh_R,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Mirollo-Strogatz');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northeast')
set(gca, 'yscale','log')
%set(gca, 'xscale', 'log')
%}

% Figure 6: Average time vs. Prob
%{
fig6 = figure(6);
clf
set(fig6, 'Position', [315, 250, 450, 380])
%yyaxis right % Times to Sync
hold on
plot(1./cycle_avg_Opt, sync_prob_Opt, 'kx', 'Linewidth', 1.0)
plot(1./cycle_avg_RL, sync_prob_RL, 'r*', 'Linewidth', 1.0)
plot(1./cycle_avg_MS_alt, sync_prob_MS_alt, 'b+', 'Linewidth', 1.0)
hold off
grid on
%axis([0.0,0.2, 4, 1000])
xlabelh = xlabel('Inverse of Average Time to \epsilon-synchronization');
ylabelh_R = ylabel('Probability of Synchronization');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh_R,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Mirollo-Strogatz');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'southeast')
%set(gca, 'yscale','log')
%set(gca, 'xscale', 'log')
%}

%% Functions
% Standard PRF for PCO Synchronization
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
% Returns:
%   y = float; the value of the PRF at phase
function [y] = StandardPRF(phase, threshold, refractory, alpha)
if phase < threshold && phase > refractory
    if phase < 0.5 * threshold
        y = -phase;
    else
        y = threshold - phase;
    end
else
    y = 0.0;
end
y = alpha * y;
end

% Second PRF for PCO Synchronization (non-identical frequencies)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
% Returns:
%   y = float; the value of the PRF at phase
function [y] = StandardPRF2(phase, threshold, refractory, alpha)
if phase < threshold && phase > refractory
    if phase < 0.5 * threshold
        y = -sqrt(2/3)*0.5*threshold*sin((pi/threshold)*phase);
    else
        y = sqrt(2/3)*0.5*threshold*sin((pi/threshold)*phase);
    end
else
    y = 0.0;
end
y = alpha * y;
end

% Third PRF for PCO Synchronization (based on RL training)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
% Returns:
%   y = float; the value of the PRF at phase
function [y] = GeneralPRF(phase, threshold, refractory, alpha)
c1 = 0.07 + 0.43*(exp(-0.23*((5)-1))); % [0.07, 0.5]
c2 = 0.77 - 0.27*(exp(-0.13*((5)-1))); % [0.5, 0.77]
if phase < threshold && phase > refractory
    if phase < c1 * threshold
        y = -phase;
    elseif phase < c2 * threshold
        y = -(c1/(1-c1))*(threshold - phase);
    else
        y = threshold - phase;
    end
else
    y = 0.0;
end
y = alpha * y;
end

% Peskin/Mirollo-Strogatz PRF for PCO Synchronization
% !! Simulation does not incorporate absorbtion
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
% Returns:
%   y = float; the value of the PRF at phase
function [y] = MSPRF(phase, threshold, refractory, alpha)
%epsilon = 0.1008;
%beta = 4;
global epsilon beta;
c = exp(beta)-1;%eeb = exp(beta*epsilon/threshold);
if phase < threshold && phase > refractory
    x = threshold * log(1+c*phase/threshold)/beta;
    phase_plus = threshold*(exp(beta*(x+epsilon)/threshold)-1)/c;
    if phase_plus > threshold
        phase_plus = threshold;
    end
    y = (phase_plus - phase);
else
    y = 0.0;
end
end

% Compute containing arc, lambda, of the PCO network
% Parameters:
%   theta = float; vector of oscillator phases
%   threshold = float; oscillator threshold
% Returns:
%   lambda = float; containing arc of theta
function [lambda] = ContainingArc(theta, threshold)
constant = ones(size(theta))*theta(1);
if isequal(theta, constant) % If the phases are all equal
%if theta == zeros(size(theta)) | theta == threshold*ones(size(theta))
%if theta <= (1E-15)*ones(size(theta)) | theta >= (threshold-(1E-15))*ones(size(theta))
    lambda = 0;
else
    buff = sort(theta); % Modulo may be unnecessary
    v = mod(buff - circshift(buff,1), threshold);
    lambda = threshold - max(v);
end
end

function [t, s, c] = PCO_Sync_Sim(PRF, theta, w_0, omega, cycles, A, alpha, D, threshold, pulse_delay, eps_sync)
k = 1; % Simulation Iteration Index
N = length(theta); % Number of oscillators in the network
% Store initial state
t = zeros(1,1); % Vector of times when data is stored
s = zeros(1,N); % Reset state values
s(k,:) = theta; % State of network at time(k)
c = zeros(1,1); % Reset arc values
c(k) = ContainingArc(s(k,:), threshold); % Containing Arc of network at time(k)

delay = cell(1,N); % Empty cell array for delay timer lists

while t(k) < cycles * (threshold/w_0)
    % First, find smallest time to next received pulse
    receivetime = threshold/(min(omega)) + eps; % Largest possible time between firings
    index_d = 0;
    j = 0;
    for i = 1:1:length(delay)
        if isempty(delay{i}) == false
            if min(delay{i}) < receivetime
                [receivetime, j] = min(delay{i});
                index_d = i;
            end
        end
    end
    % Next, find smallest time to next pulse firing
    [firingtime, index_f] = min((threshold - theta)./omega);
    % Find overall smallest time between these two values
    [timestep, index] = min([receivetime, firingtime]);
    % Updated phases to next event
    theta = theta + timestep*omega;
    % Decrease delays by timestep
    for i = 1:1:length(delay)
        if isempty(delay{i}) == false % If oscillator i has a delay timer
            delay{i} = delay{i} - timestep;
        end
    end
    % Store data
    k = k + 1;
    t(k,:) = t(k-1,:) + timestep;
    s(k,:) = theta;
    c(k) = ContainingArc(s(k,:), threshold);
    
    % Calculate change to network state based on the current event
    d_theta = zeros(size(theta));
    if index == 1 % Next event is received pulse by oscillator index_d
        % Oscillator index_d changes its phase according to the PRF
        d_theta(index_d) = PRF(theta(index_d), threshold, D, alpha);
        % Remove timer
        delay{index_d}(j) = [];
    elseif index == 2 % Next event is pulse firing by oscillator index_f
        % Oscillator index_f resets is phase.
        d_theta(index_f) = -threshold;
        for i = 1:1:N
            if A(index_f,i) == 1 % If oscillator i is connected to oscillator index_f...
                % Add delay timer to oscillator i for next received pulse
                delay{i} = [delay{i}, pulse_delay(index_f,i)]; % Fixed delay
                %delay{i} = [delay{i}, 0.02 + (0.02)*rand(1)]; % Random delay
            end
        end
    end
    
    % Updated phases based on previous event
    theta = mod(theta + d_theta, threshold);
        % Need modulo in case phase value crosses the threshold due to phase change
    % Store data
    k = k + 1;
    t(k,:) = t(k-1,:);
    s(k,:) = theta;
    c(k) = ContainingArc(s(k,:), threshold);
    
    if c(k) < eps_sync
        break % End loop; network has sufficiently synchronized
    end
end % End while loop
end % End function

% Gets the adjacency matrix for the desired network topology
% Pamameters:
%   N = int; Number of nodes/oscillators in the network
%   option = string; Desired topology type ('all', 'ring', 'dring',
%   'chain', 'dchain', 'lop2', 'lop3', 'near2', 'near3', 'far2', 'far3',
%   'near1far1', 'star', 'custom')
%   vector = int vector; List of nodes/oscillators that connect to
%   node/oscillator 1, and then shifted for remaining nodes - used with
%   'custom' option
%   symmetric = boolean; Ensures that the topology is symmetric (A = A')
% Returns:
%   A = matrix; N by N unweighted adjacency matrix
function [A] = GetTopology(N,option,vector,symmetric)
arguments % Default values for all arguments
    N = 2
    option = 'all'
    vector = []
    symmetric = false
end
switch option
    case 'all' % All-to-All
        A = MakeTopology(N, 2:N);
    case 'ring' % Bidirectional Ring
        A = MakeTopology(N, 2, true);
    case 'dring' % Directional Ring
        A = MakeTopology(N, 2);
    case 'chain' % Bidirectional Chain
        A = MakeTopology(N, 2, true); % Similar to Directional Ring
        A(N,1) = 0; A(1,N) = 0; % Include to make Chain
    case 'dchain' % Directional Chain
        A = MakeTopology(N, 2); % Similar to Directional Ring
        A(N,1) = 0; % Include to make Chain
    case 'lop2' % Lopsided, two on left
        A = MakeTopology(N, [2,3]);
    case 'lop3' % Lopsided, three on left
        A = MakeTopology(N, [2,3,4]);
    case 'near4' % Nearest two on both sides
        A = MakeTopology(N, [2,3], true);
    case 'near6' % Nearest three on both sides
        A = MakeTopology(N, [2,3,4], true);
    case 'far2' % Two farthest
        mid = round(N/2);
        A = MakeTopology(N, [mid,mid+1]);
    case 'far3' % Three farthest
        mid = round(N/2);
        A = MakeTopology(N, [mid,mid+1,mid+2]);
    case 'near2far1' % Two nearest and farthest
        mid = round(N/2);
        A = MakeTopology(N, [2,mid+1,N]);
    case 'star'
        A = MakeTopology(N, 2:N);
        A(2:N,2:N) = zeros(N-1);
    case 'custom'
        A = MakeTopology(N, vector, symmetric); % Custom
        % Include other modifications to A here.
    otherwise % Custom
        A = MakeTopology(N, [], false); % Empty
end
end

% Generates the adjacency matrix using which nodes connect to
% node/oscillator 1 and using circshift.
% Parameters:
%   N = int; Number of nodes/oscillators in the network
%   edges = int vector; List of nodes/oscillators that connect to
%   node/oscillator 1, and then shifted for remaining nodes
%   sym = boolean; Ensures that the topology is symmetric (A = A')
% Returns:
%   A = matrix; N by N unweighted adjacency matrix
function [A] = MakeTopology(N,edges,sym)
arguments % Default values for all arguments
    N = 2;
    edges (1,:) double = [];
    sym logical = false;
end
diagonal = eye(N);
A = zeros(N); % Initialize adjacency matrix
for i = edges % Add all the edges
    if i == 1
        continue % Don't make a node adjacent to itself
    elseif i > N
        break % Stop adding connections after the last node
    end
    A = A + circshift(diagonal,1-i);
end
if sym
    A = double((A + A.')>0); % Make the graph symmetric
end
end

function DisplayPercentageComplete(iter,total,time,step)
percent = 100*iter/total;
for i = 0:step:100
    if percent == i
        disp([num2str(i),'% Complete in ',num2str(time),' seconds'])
        break % We done
    end
end
end

% Check that the given adjacency matrix defines a connected network
% Parameters:
%   A = matrix; square unweighted adjacency matrix - Use MakeTopology or
%   GetTopology to get this matrix
function CheckTopology(A)
% Construct graph and check if connected.
G=digraph(A);
bins_weak = conncomp(G,'Type','weak');
bins_strong = conncomp(G,'Type','strong');
if all(bins_strong == 1)
    disp("Graph is strongly connected.")
elseif all(bins_weak == 1)
    disp("Graph is weakly connected.")
else % Something went wrong
    disp("Graph is not connected. Press ENTER to continue.")
    pause
end
end

function [score] = Score_Population(time, arcs, weights, epsilon)
total1 = 0;
total2 = 0;
for i = 1:1:(0.5*(length(time)-1))
    delta_time = time(2*i) - time(2*i - 1);
    C_avg = 0.5*(arcs(2*i-1) + arcs(2*i));
    total1 = total1 + (C_avg*delta_time);
    % Potentially incorporate additional scoring criteria
    if arcs((2*i)-1) < epsilon && arcs(2*i) < epsilon
        % Arcs is under the threshold for the entire delta_time
        total2 = total2 + delta_time;
    elseif arcs((2*i)-1) < epsilon || arcs(2*i) < epsilon
        % Arcs is under the threshold some fraction of delta_time
        arcmin = min(arcs((2*i)-1:(2*i)));
        arcmax = max(arcs((2*i)-1:(2*i)));
        frac = (epsilon - arcmin)/(arcmax-arcmin);
        total2 = total2 + frac*delta_time;
    %else Arcs is not under the threshold for the entire delta_time
    end
end
score = weights(1)/total1 + weights(2)*total2;
end
