% PCO Simulation
% Based on PCO_Sim060419
% Calculates estimated containing arc for each oscillator
% Last Modified: 7/6/2020
% Issues with state estimates when estimates get very close to synchronization

%% Simulation Parameters
cycles = 250; % Number of evolution cycles
alpha = 0.1; % PCO coupling strength - (0,1]
threshold = 2*pi; % Oscillator phase threshold
D = 0.0; % Refractory Period
N = 11; % Number of PCOs (integer >= 2)
A_option = 'near4';
A = GetTopology(N,A_option);
%A = GetTopology(N,'custom',[2:5],true); % Custom topology
%rng('default');
w_0 = 2*pi; % Nominal oscillator evolution frequency
omega = w_0 * ones(1,N) + (0.01)*2*(rand(1,N)-0.5); % vector of oscillator frequencies
rand_delay = (0.02)*ones(N,N) + (0.08)*rand(N,N); % Random matrix
pulse_delay = 0.5 * (rand_delay + rand_delay'); % (symmetric) matrix of pulse propogation delays
% Nominal (simulation) time to complete simulation (in seconds): cycles * (threshold/w_0);
theta_int = 1;
switch theta_int
    case 1 % Random
        %theta = sort(rand(1,N),'descend')*threshold; % vector of initial node phases [0,1]
        theta = rand(1,N)*threshold; % vector of (unsorted) initial node phases [0,1]
    case 2 % largest possible spacing
        spacing = threshold/(1*N);
        theta = ((N-1)*spacing:-spacing:0) + 0.01;
    case 3 % Use last initial state
        theta = theta_start;
    otherwise % Manually determined
        %theta = 2*pi*[0.51, 0.45, 0.41, 0.33, 0.21, 0.03];
        
        % Examples of Standard PRF not synchronizing at alpha \approx 0.1
        %theta = [0.4061, 0.4463, 2.4789, 2.7406, 5.1939, 5.5778];
        %theta = [0.4947, 1.3557, 2.3719, 3.7880, 4.1779, 5.8626];
        %theta = [0.2614, 0.5858, 1.6658, 2.0519, 4.3616, 4.7707];
        %theta = [0.4789, 0.6781, 2.6675, 3.0400, 3.5836, 6.1088];
        %theta = [0.9281, 1.0085, 2.4833, 3.2229, 4.6989, 5.6939];
        %theta = [1.3712, 6.076, 2.7268, 4.9312, 3.3001, 2.0815];
        theta = [1.2275, 0.5185, 4.5723, 3.2763, 3.2467, 0.1608];
end
% Network state validation
if length(theta) ~= N
    disp("Error: Theta Size Mismatch")
    pause;
end
theta_start = theta; % Store initial oscillator states (for additional simulations)

%% Simulation
% Run simulation and return results
[time, state, C_Arc, C_Arc_est] = PCO_Sync_Sim(@StandardPRF, theta, w_0, omega, cycles, A,...
                                                alpha, D, threshold, pulse_delay);
[time2, state2, C_Arc2, C_Arc_est2] = PCO_Sync_Sim(@MSPRF, theta, w_0, omega, cycles, A,...
                                                alpha, D, threshold, pulse_delay);

%% Plots
% Figure 1: PCO phases over time
%{
fig1 = figure(1);
clf
set(fig1, 'Position', [10, 300, 500, 340])
hold on
%plot(time,state, 'Linewidth', 1.5)
plot(time2,state2, 'Linewidth', 1.5)
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('PCO Phases (\theta)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,cycles*(threshold/w_0),0,threshold])
%}

% Figure 2: Containing Arc over time
%
fig2 = figure(2);
clf
set(fig2, 'Position', [615, 400, 500, 240])
hold on
for i = 1:N
    plot(time,C_Arc_est(:,i), 'Linewidth', 1.5)
end
plot(time,C_Arc, 'k-', 'Linewidth', 1.5)
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('Containing Arc (\Lambda)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,cycles*(threshold/w_0),0,threshold*(N-1)/N])
%set(gca, 'yscale','log')
%}

% Figure 2: Containing Arc over time
%
fig3 = figure(3);
clf
set(fig3, 'Position', [615, 50, 500, 240])
hold on
for i = 1:N
    plot(time2,C_Arc_est2(:,i), 'Linewidth', 1.5)
end
plot(time2,C_Arc2, 'k-', 'Linewidth', 1.5)
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('Containing Arc (\Lambda)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,cycles*(threshold/w_0),0,threshold*(N-1)/N])
%set(gca, 'yscale','log')
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
end

% Third PRF for PCO Synchronization (based on RL training)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
% Returns:
%   y = float; the value of the PRF at phase
function [y] = GeneralPRF(phase, threshold, refractory, alpha)
c1 = 0.33; % [0.14, 0.5]
c2 = 0.62; % [0.5, 0.8]
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
epsilon = 0.02;
beta = 5;
c = exp(beta)-1;
%eeb = exp(beta*epsilon/threshold);
if phase < threshold && phase > refractory
    x = threshold * log(1+c*phase/threshold)/beta;
    phase_plus = threshold*(exp(beta*(x+epsilon)/threshold)-1)/c;
    %phase_plus = phase*eeb + (threshold/c)*(eeb - 1);
    if phase_plus > threshold
        phase_plus = threshold;
    end
    y = (phase_plus - phase)/alpha;
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

function [t, s, c, c_est] = PCO_Sync_Sim(PRF, theta, w_0, omega, cycles, A, alpha, D, threshold, pulse_delay)
k = 1; % Simulation Iteration Index
N = length(theta); % Number of oscillators in the network
% Store initial state
t = zeros(1,1); % Vector of times when data is stored
s = zeros(1,N); % Reset state values
s(k,:) = theta; % State of network at time(k)
c = zeros(1,1); % Reset arc values
c(k) = ContainingArc(s(k,:), threshold); % Containing Arc of network at time(k)
delta_matrix = cell(1,N); % Current estimate of phase differences from perspective of each oscillator
    % delta_theta_est{i} is oscillator i's estimate of the phase
    % differences to known oscillators from the previous cycle
m = zeros(1,N);
c_est = zeros(1,N); % Containing arc of estimated network state for oscillators
    % c_est(k,i) is the estimate of the containing arc from oscillator i's perspective at timestep k
[theta_sort, osc_order] = sort(theta,'descend');
for i = 1:N
    % Assume perfect initial knowledge of state
    delta_matrix{osc_order(i)} = zeros(1,sum(A(osc_order(i),:)));
    count = 0;
    for j = circshift(osc_order,-i)
        if A(osc_order(i),j) == 1
            count = count + 1;
            delta_matrix{osc_order(i)}(count) = mod(theta(j) - theta_sort(i), threshold);
            if theta(j) < theta_sort(i)
                m(osc_order(i)) = m(osc_order(i)) + 1; % Initial value for the number of oscillator pulses received after own firing
            end
        end
    end
    c_est(k,osc_order(i)) = ContainingArc([delta_matrix{i}, 0],threshold);
end
% Display stuff for debugging
%{
disp(m)
for i = 1:N
    disp(delta_matrix{i})
end
pause
%}
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
    for i = 1:N
        c_est(k,i) = ContainingArc([delta_matrix{i}, 0],threshold);
    end
    
    % Calculate change to network state based on the current event
    d_theta = zeros(size(theta));
    if index == 1 % Next event is received pulse by oscillator index_d
        % Oscillator index_d changes its phase according to the PRF
        d_theta(index_d) = alpha * PRF(theta(index_d), threshold, D, alpha);
        % Remove timer
        delay{index_d}(j) = [];
        
        % Oscillator index_d updates phase difference estimates of other oscillators
        % m is the number of received pulses since its own firing
        m(index_d) = m(index_d) + 1; % Increment received pulse counter
        delta_matrix{index_d}(m(index_d)) = mod(-(theta(index_d) + d_theta(index_d)),threshold); % Update for oscillator that fired
        for i = 1:length(delta_matrix{index_d})
            if i == m(index_d)
                %delta_matrix{index_d}(i) = mod(-(theta(index_d) + d_theta(index_d)),threshold); % Update for oscillator that fired
            else %if i ~= m(index_d) % For all other oscillators; assume pulse was received at the same time
                % Change to phase difference estimate
                stuff = alpha * PRF(mod(theta(index_d)+delta_matrix{index_d}(i), threshold), threshold, D, alpha) - d_theta(index_d);
                delta_matrix{index_d}(i) = delta_matrix{index_d}(i) + stuff;
            end
        end
        
    elseif index == 2 % Next event is pulse firing by oscillator index_f
        % Oscillator index_f resets is phase.
        d_theta(index_f) = -threshold;
        for i = 1:1:N
            if A(index_f,i) == 1 % If oscillator i is connected to oscillator index_f...
                % Add delay timer to oscillator i for next received pulse
                %delay{i} = [delay{i}, pulse_delay(index_f,i)]; % Fixed delay
                delay{i} = [delay{i}, 0.02 + (0.08)*rand(1)]; % Random delay
            end
        end
        
        if m(index_f) < length(delta_matrix{index_f}) % If pulses received is less than state differences that oscilator index_f is tracking...
            % Remove excess state differences
            delta_matrix{index_f}(m(index_f)+1:end) = [];
        end % Otherwise, leave delta_matrix alone
        m(index_f) = 0; % reset received pulse counter
        
        % Oscillator index_f updates phase difference estimates of all oscillators
        for i = 1:length(delta_matrix{index_f})
            % Change to phase difference estimate
            stuff = alpha * PRF(mod(delta_matrix{index_f}(i), threshold), threshold, D, alpha);
            delta_matrix{index_f}(i) = delta_matrix{index_f}(i) + stuff;
        end
        
    end
    
    % Updated phases based on previous event
    theta = mod(theta + d_theta, threshold);
        % Need modulo in case phase value crosses the threshold due to phase change
    %theta_est = theta_est + d_theta_est; % Might need modulo...
    % Store data
    k = k + 1;
    t(k,:) = t(k-1,:);
    s(k,:) = theta;
    c(k) = ContainingArc(s(k,:), threshold);
    for i = 1:N
        c_est(k,i) = ContainingArc([0, delta_matrix{i}],threshold);
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
