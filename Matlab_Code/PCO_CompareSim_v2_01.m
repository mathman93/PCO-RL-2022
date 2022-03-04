% PCO Simulation Comparison - Version 2.01
% Based on Version 2
% Compares performance of PCO synchronization PRFs/algorithms
% Focuses on generating specific comparison graph on amount of synchronization in non-ideal enviroments.
% Last Modified: 12/1/2021

%% Simulation Parameters
cycles = 151; % Maximum number of evolution cycles
N = 10; % Number of PCOs (integer >= 2)
%alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999];
alphas = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999];
%alphas = [0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40];
threshold = 2*pi; % Oscillator phase threshold
w_0 = 2*pi; % Nominal oscillator evolution frequency
D = 0.0; % Refractory Period
A_option = 'ERGrand';
%A = GetTopology(N,A_option);
%A = GetTopology(N,'custom',[2:5],true); % Custom topology
%A = GetTopology(N,'ERGrand',0.3,true); % Random topology
out_degree = sum(A,2); % Outdegree of oscillators
out_degree_avg = sum(out_degree)/N; % Average network outdegree
in_degree = sum(A,1) % Indegree of oscillators
in_degree_avg = sum(out_degree)/N; % Average network indegree
eps_sync = 0.001*threshold; % epsilon-synchronization threshold

%theta_start = theta; % Store initial oscillator states (for additional simulations)

area_1 = (0.5*threshold)^2;
c = [0.095 + 0.405*(exp(-0.29*((in_degree')-1))), 0.755 - 0.255*(exp(-0.14*((in_degree')-1)))]; % learned PRF parameters
area_2 = 0.5*threshold^2 * (c(:,1).^2 + (1-c(:,2)).^2 + (1-c(:,1)).*c(:,1) - c(:,1).*(1-c(:,2)).*(1-c(:,2))./(1-c(:,1)));

% MS parameters
epsilon = 0.05*ones(N,1);
beta = 5*ones(N,1);
% Estimate value of epsilon that gives same "area" as GeneralPRF
%{
p_step = threshold/1000;
p = 0:p_step:threshold;
stateMS = log(1+(exp(beta)-1)*(p./threshold))./beta;
for i=1:10
    phaseMS = threshold.*(exp(beta.*(stateMS + epsilon))-1)./(exp(beta)-1);
    p_plus = min((phaseMS - p),(threshold - p));
    area_3 = p_step * (sum(p_plus,2) - 0.5*(p_plus(:,1)+p_plus(:,end)));
    epsilon = epsilon.*(area_2.*(1.*alphas(end))./area_3);
end
disp(epsilon)
%}

%% Simulation
rng('default');
sample_size = 2000; % Number of simulations to run/compare
w_var = 0.05; % oscillator frequency variation (percentage of w_0)
tau_min = 0.01; % minimum pulse delay (fraction of cycle)
tau_max = 0.08; % maximum pulse delay (fraction of cycle)
r_max = 0.5; % partition parameter for Javed2021/Nishimura2012
r_min = 0.5; % partition parameter for Javed2021/Nishimura2012
r_part = r_min*ones(N,1) + (r_max - r_min)*rand(N,1);
L = zeros(sample_size,1); % Initial containing arc of oscillator for each set of simulations
maxTime = zeros(sample_size,4); % Time to synchronize for each simulation
SyncAmount = zeros(sample_size,4); % Average containing arc value once network reaches below eps_sync.
sample = 0; % Sample number index
thetas = zeros(sample_size,N);
omegas = zeros(sample_size,N);

while sample < sample_size % Generate initial phase values for simulations
    i = sample + 1;
    thetas(i,:) = rand(1,N)*threshold; % vector of (unsorted) initial node phases [0,1]
    % Run simulation and return results
    Lambda_0 = ContainingArc(thetas(i,:), threshold); % Initial containing arc
    %
    if Lambda_0 < 0.5*threshold
        continue % Don't include results for "small" initial thresholds
    end
    %}
    omegas(i,:) = w_0 * (ones(1,N) + (w_var)*2*(rand(1,N)-0.5)); % vector of oscillator frequencies
    L(i) = Lambda_0;
    sample = sample + 1; % Increment sample counter
end

C_avg_Opt = zeros(sample_size, length(alphas));
C_avg_RL = zeros(sample_size, length(alphas));
C_avg_Kling = zeros(sample_size, length(alphas));
C_avg_Nishi = zeros(sample_size, length(alphas));
C_avg_MS = zeros(sample_size, length(alphas));
for ai = 1:length(alphas)
alpha = alphas(ai)*ones(N,1);
disp('Coupling Strength = ' + string(alphas(ai)))
start_time = tic; % For timing the simulations
for sample = 1:sample_size
    theta = thetas(sample,:);
    omega = omegas(sample,:);
    %
    [time1, state1, C_Arc1, t_sync1, C_avg1] = PCO_Sync_Sim(@StandardPRF, theta, w_0, omega, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    [time2, state2, C_Arc2, t_sync2, C_avg2] = PCO_Sync_Sim(@GeneralPRF, theta, w_0, omega, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    %
    [time3, state3, C_Arc3, t_sync3, C_avg3] = PCO_Sync_Sim(@KlinglmayrPRF, theta, w_0, omega, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    %
    [time4, state4, C_Arc4, t_sync4, C_avg4] = PCO_Sync_Sim(@NishimuraPRF, theta, w_0, omega, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    %
    [time5, state5, C_Arc5, t_sync5, C_avg5] = PCO_Sync_Sim(@MSPRF, theta, w_0, omega, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    % Store values for analysis later
    C_avg_Opt(sample, ai) = C_avg1;
    C_avg_RL(sample, ai) = C_avg2;
    C_avg_Kling(sample, ai) = C_avg3;
    C_avg_Nishi(sample, ai) = C_avg4;
    C_avg_MS(sample, ai) = C_avg5;
    % Calculate percentage complete and display
    DisplayPercentageComplete(sample, sample_size, toc(start_time), 1);
end % for sample
end % for ai

%% Plots
black = [0,0,0]; red = [1,0,0]; blue = [0,0,1]; purp = [0.6,0,1]; green = [0,0.7,0.2]; % color defs
% Figure 1: PCO phases over time
%{
fig1 = figure(1);
clf
set(fig1, 'Position', [10, 300, 500, 340])
hold on
%plot(time1,state1, 'Linewidth', 1.5)
plot(time4,state4, 'Linewidth', 1.5)
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
plot(time1,C_Arc1, '-', 'Color', black, 'Linewidth', 1.5) % StandardPRF
plot(time2,C_Arc2, '-', 'Color', red, 'Linewidth', 1.5) % GeneralPRF
plot(time3,C_Arc3, '-', 'Color', blue, 'Linewidth', 1.5) % Klinglmayr2017
plot(time4,C_Arc4, '-', 'Color', purp, 'Linewidth', 1.5) % Nishimura2012
plot(time5,C_Arc5, '-', 'Color', green, 'Linewidth', 1.5) % Mirollo-Strogatz1990
%plot(time1,C_avg1*ones(size(time1)), '--', 'Color', black, 'Linewidth', 1.5)
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('Containing Arc (\Lambda)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Wang-Doyle (2012)', 'Learned PRF', 'Klinglmayr et.al. (2017)', 'Nishimura-Friedman (2011)', 'Mirollo-Strogatz (1990)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northeast')
%axis([0,max(maxTime(end,:)),0,threshold*(N-1)/N])
set(gca, 'yscale','log')
%}

%% PRF plots
p_steps = 1000;
p = linspace(0,threshold,p_steps+1);
StandF = zeros(size(p));
RLF = zeros(size(p));
KlingF = zeros(size(p));
%JavedF = zeros(N,length(p));
NishiF = zeros(size(p));
MSF = zeros(size(p));
r_part = r_min*ones(N,1) + (r_max - r_min)*rand(N,1);
alpha = 0.05*ones(N,1);
params = [alpha, epsilon, beta, tau_min*ones(N,1), tau_max*ones(N,1), c(:,1), c(:,2), r_part];
for i = 1:length(p)
    StandF(i) = StandardPRF(p(i),1,threshold,D,params);
    RLF(i) = GeneralPRF(p(i),1,threshold,D,params);
    KlingF(i) = KlinglmayrPRF(p(i),1,threshold,D,params);
    NishiF(i) = NishimuraPRF(p(i),1,threshold,D,params);
    MSF(i) = MSPRF(p(i),1,threshold,D,params);
end
KlingF_mod = mod(KlingF + (threshold/2), threshold) - (threshold/2);
KlingF_area = (threshold/p_steps) * (sum(abs(KlingF_mod),2) - 0.5*(abs(KlingF_mod(:,1))+abs(KlingF_mod(:,end))));
%disp(KlingF_area/area_1)
NishiF_area = (threshold/p_steps) * (sum(abs(NishiF),2) - 0.5*(abs(NishiF(:,1))+abs(NishiF(:,end))));
%disp(NishiF_area/area_1)
MSF_area = (threshold/p_steps) * (sum(abs(MSF),2) - 0.5*(abs(MSF(:,1))+abs(MSF(:,end))));
%disp(MSF_area/area_1)

% Figure 9: PRF comparison
%{
fig9 = figure(9);
clf
set(fig9, 'Position', [55, 50, 350, 440])
hold on
plot(p, StandF, 'k-', 'LineWidth', 1)
plot(p, RLF, 'r--', 'LineWidth', 1)
plot(p, KlingF, 'b-', 'LineWidth', 1)
%plot(p, KlingF_mod, 'c-', 'LineWidth', 1)
plot(p, NishiF, '-.', 'LineWidth', 1)
plot(p, MSF, '-', 'Linewidth', 1)
hold off
grid on
xlabelh = xlabel('Oscillator Phase');
ylabelh = ylabel('Phase Response Function');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF (l=1.0)', 'Learned PRF (l=1.0)', 'Klinglmayr (2017)', 'Nishimura (2011)', 'Mirollo-Strogatz (1990)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'southwest')
axis([0,threshold, -threshold,threshold/2])
xticks([0 (pi/2) pi (3*pi/2) 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
yticks([-2*pi -3*pi/2 -pi (-pi/2) 0 (pi/2) pi 3*pi/2 2*pi])
yticklabels({'-2\pi', '-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
%set(gca, 'yscale','log')
%}

%% Comparison Plots
% Check if test synchronized
sync_bool_Opt = C_avg_Opt < 0.5*threshold;
sync_bool_RL = C_avg_RL < 0.5*threshold;
sync_bool_K = C_avg_Kling < 0.5*threshold;
sync_bool_N = C_avg_Nishi < 0.5*threshold;
sync_bool_MS = C_avg_MS < 0.5*threshold;
% Pre-allocate variables
med_Opt = zeros(size(alphas));
mean_Opt = zeros(size(alphas));
max_Opt = zeros(size(alphas));
min_Opt = zeros(size(alphas));
std_Opt = zeros(size(alphas));
sync_prob_Opt = zeros(size(alphas));
med_RL = zeros(size(alphas));
mean_RL = zeros(size(alphas));
max_RL = zeros(size(alphas));
min_RL = zeros(size(alphas));
std_RL = zeros(size(alphas));
sync_prob_RL = zeros(size(alphas));
med_K = zeros(size(alphas));
mean_K = zeros(size(alphas));
max_K = zeros(size(alphas));
min_K = zeros(size(alphas));
std_K = zeros(size(alphas));
sync_prob_K = zeros(size(alphas));
med_N = zeros(size(alphas));
mean_N = zeros(size(alphas));
max_N = zeros(size(alphas));
min_N = zeros(size(alphas));
std_N = zeros(size(alphas));
sync_prob_N = zeros(size(alphas));
med_MS = zeros(size(alphas));
mean_MS = zeros(size(alphas));
max_MS = zeros(size(alphas));
min_MS = zeros(size(alphas));
std_MS = zeros(size(alphas));
sync_prob_MS = zeros(size(alphas));
for ai = 1:length(alphas)
    % Standard PRF
    %data = C_avg_Opt(sync_bool_Opt(:,ai), ai);
    data = C_avg_Opt(:, ai);
    med_Opt(ai) = median(data);
    mean_Opt(ai) = mean(data);
    max_Opt(ai) = max(data);
    min_Opt(ai) = min(data);
    std_Opt(ai) = std(data);
    sync_prob_Opt(ai) = sum(sync_bool_Opt(:,ai))/sample_size;
    % RL PRF
    %data = C_avg_RL(sync_bool_RL(:,ai), ai);
    data = C_avg_RL(:, ai);
    med_RL(ai) = median(data);
    mean_RL(ai) = mean(data);
    max_RL(ai) = max(data);
    min_RL(ai) = min(data);
    std_RL(ai) = std(data);
    sync_prob_RL(ai) = sum(sync_bool_RL(:,ai))/sample_size;
    % Klinglmayr PRF
    %data = C_avg_Kling(sync_bool_K(:,ai), ai);
    data = C_avg_Kling(:, ai);
    med_K(ai) = median(data);
    mean_K(ai) = mean(data);
    max_K(ai) = max(data);
    min_K(ai) = min(data);
    std_K(ai) = std(data);
    sync_prob_K(ai) = sum(sync_bool_K(:,ai))/sample_size;
    % Nishimura PRF
    %data = C_avg_Nishi(sync_bool_N(:,ai), ai);
    data = C_avg_Nishi(:, ai);
    med_N(ai) = median(data);
    mean_N(ai) = mean(data);
    max_N(ai) = max(data);
    min_N(ai) = min(data);
    std_N(ai) = std(data);
    sync_prob_N(ai) = sum(sync_bool_N(:,ai))/sample_size;
    % Mirollo-Strogatz PRF
    %data = C_avg_MS(sync_bool_MS(:,ai), ai);
    data = C_avg_MS(:, ai);
    med_MS(ai) = median(data);
    mean_MS(ai) = mean(data);
    max_MS(ai) = max(data);
    min_MS(ai) = min(data);
    std_MS(ai) = std(data);
    sync_prob_MS(ai) = sum(sync_bool_MS(:,ai))/sample_size;
end

% Figure 3: Compare Sync steady-state
%{
fig3 = figure(3);
clf
set(fig3, 'Position', [55, 50, 440, 400])
hold on
errorbar(alphas, mean_Opt, std_Opt, 'kx--', 'LineWidth', 1.4)
errorbar(alphas, mean_RL, std_RL, 'rx--', 'LineWidth', 1.4)
errorbar(KlingF_area/area_1, mean_K, std_K, 'bx', 'LineWidth', 1.4)
errorbar(NishiF_area/area_1, mean_N, std_N, 'x', 'LineWidth', 1.4)
errorbar(MSF_area/area_1, mean_MS, std_MS, 'x', 'LineWidth', 1.4)
hold off
grid on
xlabelh = xlabel('Coupling Strength, l');
ylabelh = ylabel('Steady-state Containing Arc');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Klinglmayr et.al. (2017)', 'Nishimura (2011)', 'Mirollo-Strogatz (1990)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northeast')
%axis([0,1.0, 0.03,0.26])
%set(gca, 'yscale','log')
%}

% Figure 3-Alt: Compare Sync steady-state
%
fig3 = figure(3);
clf
set(fig3, 'Position', [55, 50, 440, 320])
hold on
%{
errorbar(alphas, mean_Opt, mean_Opt-min_Opt, max_Opt-mean_Opt, 'x--', 'Color', black, 'LineWidth', 1.4)
errorbar(alphas, mean_RL, mean_RL-min_RL, max_RL-mean_RL, '*--', 'Color', red, 'LineWidth', 1.4)
errorbar(alphas, mean_K, mean_K-min_K, max_K-mean_K, '^--', 'Color', blue, 'LineWidth', 1.4)
errorbar(alphas, mean_N, mean_N-min_N, max_N-mean_N, 'v--', 'Color', purp, 'LineWidth', 1.4)
errorbar(alphas, mean_MS, mean_MS-min_MS, max_MS-mean_MS, 'o--', 'Color', green, 'LineWidth', 1.4)
%}
%
plot(alphas, mean_RL, '*--', 'Color', red, 'LineWidth', 1.4)
plot(alphas, mean_K, '^--', 'Color', blue, 'LineWidth', 1.4)
plot(alphas, mean_Opt, 'x--', 'Color', black, 'LineWidth', 1.4)
plot(alphas, mean_N, 'v--', 'Color', purp, 'LineWidth', 1.4)
plot(alphas, mean_MS, 'o--', 'Color', green, 'LineWidth', 1.4)
%}
hold off
grid on
xlabelh = xlabel('Coupling Strength, $l$', 'interpreter', 'latex');
ylabelh = ylabel('Steady-state Containing Arc', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Learned PRF', 'Klinglmayr et.al. (2017)', 'Wang-Doyle (2012)', 'Nishimura-Friedman (2011)', 'Mirollo-Strogatz (1990)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',11);%, 'Position', [0.43, 0.43, 0.45, 0.25])
axis([0,1.0, 0.0,4.5])
%set(gca, 'yscale','log')
% Reorder plots
%chi = get(gca, 'Children');
%set(gca, 'Children', circshift(chi,2));
%}

% Figure 4: Compare Prob. Sync
%
fig4 = figure(4);
clf
set(fig4, 'Position', [555, 50, 440, 320])
hold on
plot(alphas, sync_prob_RL, '*-', 'Color', red, 'LineWidth', 1.5)
plot(alphas, sync_prob_K, '^-', 'Color', blue, 'LineWidth', 1.5)
plot(alphas, sync_prob_Opt, 'x-', 'Color', black, 'LineWidth', 1.5)
plot(alphas, sync_prob_N, 'v-', 'Color', purp, 'LineWidth', 1.5)
plot(alphas, sync_prob_MS, 'o-', 'Color', green, 'LineWidth', 1.5)
hold off
grid on
xlabelh = xlabel('Coupling Strength, $l$', 'interpreter', 'latex');
ylabelh = ylabel('Probability of Synchronization', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
%legendh = legend('Standard PRF', 'Learned PRF');
legendh = legend('Learned PRF', 'Klinglmayr et.al. (2017)', 'Wang-Doyle (2012)', 'Nishimura-Friedman (2011)', 'Mirollo-Strogatz (1990)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',11, 'Location', 'east')
%axis([0.0,1.0, 0.4, 1.0])
%set(gca, 'yscale','log')
%}

%% Functions
% Standard PRF for PCO Synchronization
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = StandardPRF(phase, index, threshold, refractory, params)
if phase < threshold && phase > refractory
    if phase < 0.5 * threshold
        y = -phase;
    else
        y = threshold - phase;
    end
else
    y = 0.0;
end
y = params(index,1) * y;
end

% Second PRF for PCO Synchronization (non-identical frequencies)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = StandardPRF2(phase, index, threshold, refractory, params)
if phase < threshold && phase > refractory
    if phase < 0.5 * threshold
        y = -sqrt(2/3)*0.5*threshold*sin((pi/threshold)*phase);
    else
        y = sqrt(2/3)*0.5*threshold*sin((pi/threshold)*phase);
    end
else
    y = 0.0;
end
y = params(index, 1) * y;
end

% Third PRF for PCO Synchronization (based on RL training)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = GeneralPRF(phase, index, threshold, refractory, params)
c1 = params(index, 6);%0.07 + 0.43*(exp(-0.23*((5)-1))); % [0.07, 0.5]
c2 = params(index, 7);%0.77 - 0.27*(exp(-0.13*((5)-1))); % [0.5, 0.77]
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
y = params(index, 1) * y;
end

% Peskin/Mirollo-Strogatz PRF for PCO Synchronization
% !! Simulation does not incorporate absorbtion
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = MSPRF(phase, index, threshold, refractory, params)
epsilon = params(index, 2);
beta = params(index, 3);
%global epsilon beta;
c = exp(beta)-1;%eeb = exp(beta*epsilon/threshold);
if phase < threshold && phase > refractory
    x = threshold * log(1+c*phase/threshold)/beta;
    phase_plus = threshold*(exp(beta*(x+epsilon)/threshold)-1)/c;
    if phase_plus > threshold
        phase_plus = threshold;
    end
    %y = (phase_plus - phase);
    y = params(index,1)*(phase_plus - phase); % Include coupling strength parameter
else
    y = 0.0;
end
end

% Klinglmayr 2017 PRF for PCO Synchronization
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = KlinglmayrPRF(phase, index, threshold, refractory, params)
tau_min = params(index, 4)*threshold;
tau_max = params(index, 5)*threshold;
if phase < refractory
    y = 0.0;
else
    p = mod(phase - tau_min, threshold);
    if p > 0.5*threshold
        %H_tilde = 0.46*p + 0.54*threshold; % From Klinglmayr2017
        H_tilde = (0.75*threshold + tau_max - tau_min) + ((p - 0.5*threshold)/(0.5*threshold))*(0.25*threshold - tau_max + tau_min); % h2
    elseif p > tau_max
        %H_tilde = 0.3261*p + 0.027*threshold; % From Klinglmayr2017
        H_tilde = tau_max + ((p-tau_max)/(0.5*threshold - tau_max))*(0.25*threshold - tau_min - 2*tau_max); % h1
    else
        H_tilde = p;
    end
    H = mod(H_tilde + tau_min, threshold);
    %y = H - phase;
    % Complicated mess to include coupling strength parameter
    H_norm = mod(H-phase + (0.5*threshold), threshold) - (0.5*threshold); % normalize
    H_scale = params(index,1) * H_norm; % scale normalized function
    y = mod(H_scale + phase, threshold) - phase; % transform back
end
end

% Nishimura 2011 PRF for PCO Synchronization (strong type II)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = NishimuraPRF(phase, index, threshold, refractory, params)
r = params(index, 8); % Might not need this.
% B0 = B1 = r
tau = params(index, 5)*threshold;
kappa = 0.02*threshold;
if phase > r*threshold % Something excitory (positive)
    y = (threshold/pi)*sin((pi/threshold)*phase);
    %y = params(index,1)*(threshold/pi)*sin((pi/threshold)*phase);
    %y = params(index,1)*(threshold - phase);
elseif phase > tau + kappa % Something inhibitory (negative) enough
    y = min(-(tau+kappa), -(threshold/pi)*sin((pi/threshold)*phase));
    %y = min(-(tau+kappa), -params(index,1)*(threshold/pi)*sin((pi/threshold)*phase));
    %y = min(-(tau+kappa), -params(index,1)*phase);
elseif phase > refractory
    y = -1 * phase;
else
    y = 0.0;
end
y = params(index,1)*y; % Include coupling strength parameter
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

function [t, s, c, t_eps, c_avg] = PCO_Sync_Sim(PRF, theta, w_0, omega, cycles, A, threshold, D, ...
    p_send, p_receive, alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r, eps_sync)
N = length(theta); % Number of oscillators in the network
PRF_params = [alpha, epsilon, beta, tau_min*ones(N,1), tau_max*ones(N,1), RL_c1, RL_c2, r];
p_s = p_send; % probability of sending a pulse when resetting (Klinglmayr)
p_r = p_receive; % probability of receiving a sent pulse (Javed)

roll_window = 10;
k = 1; % Simulation Iteration Index
k_r = 1; % Rolling Average start index
k_start = 1;
% Store initial state
t = zeros(1,1); % Vector of times when data is stored
s = zeros(1,N); % Reset state values
s(k,:) = theta; % State of network at time(k)
c = zeros(1,1); % Reset arc values
c(k) = ContainingArc(s(k,:), threshold); % Containing Arc of network at time(k)
c_roll = zeros(1,1); % Rolling Average of Containing Arc "c"
c_roll(k) = c(k);

delay = cell(1,N); % Empty cell array for delay timer lists
breakloop = false; % boolean for finishing simulation loop
end_time = cycles * (threshold/w_0);
while t(k) < end_time
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
    s(k,:) = theta; % State before action
    c(k) = ContainingArc(s(k,:), threshold);
    
    % Calculate change to network state based on the current event
    d_theta = zeros(size(theta));
    if index == 1 % Next event is received pulse by oscillator index_d
        % Oscillator index_d changes its phase according to the PRF
        d_theta(index_d) = PRF(theta(index_d), index_d, threshold, D, PRF_params);
        % Remove timer
        delay{index_d}(j) = [];
    elseif index == 2 % Next event is pulse firing by oscillator index_f
        % Oscillator index_f resets is phase.
        d_theta(index_f) = -threshold;
        if rand(1) < p_s % Send pulse with probability p_send
            for i = 1:1:N % For each oscillator in the network
                if A(index_f,i) == 1 % If oscillator i is connected to oscillator index_f...
                    % Add delay timer to oscillator i for next received pulse
                    if rand(1) < p_r % Receive pulse with probability p_receive
                        delay{i} = [delay{i}, tau_min + (tau_max-tau_min)*rand(1)]; % Random delay
                    end
                end
            end
        end
    end
    
    % Updated phases based on previous event
    theta = mod(theta + d_theta, threshold);
        % Need modulo in case phase value crosses the threshold due to phase change
    % Store data
    k = k + 1;
    t(k,:) = t(k-1,:);
    s(k,:) = theta; % State after action
    c(k) = ContainingArc(s(k,:), threshold);

    % Rolling Average of Containing Arc
    while t(k,:) > (t(k_r) + roll_window*(threshold/w_0))
        k_r = k_r + 2; % Update Rolling average start index
        roll = 0; % Reset rolling average
        for id = k_r:2:k-2 % Calculate rolling average of "c" from k_r to k
            val = 0.5*(c(id)+c(id+1)) * (t(id+1)-t(id));
            roll = roll + val;
        end
        %disp(k_r)
        c_roll(k_r-1) = c_roll(k_r-2);
        c_roll(k_r) = roll/(t(k)-t(k_r));
    end
    if k_r > 1 % Don't count consider the initially small window
        while t(k_r,:) > (t(k_start) + roll_window*(threshold/w_0))
            k_start = k_start + 1;
        end
        %
        % Condition: Change in rolling average is small
        if c(k_r) < (threshold/2) && max(c_roll(k_start:k_r)) - min(c_roll(k_start:k_r)) < eps_sync && breakloop == false
            t_eps = t(k_start,:); % Time at beginning of window
            k_eps = k_r; % For containing arc value at end of window
            end_time = t(k_r,:) + (2*roll_window)*(threshold/w_0); % set simulation to go on for 10 more seconds/cycles
            breakloop = true; % Break simulation loop after those 10 seconds/cycles
        end
        %}
        %{
        % Condition: Rolling average is small
        if c_roll(k_r) < eps_sync && breakloop == false
            t_eps = t(k_start,:);
            k_eps = k_r;
            end_time = t(k_r,:) + (2*roll_window)*(threshold/w_0); % set simulation to go on for 10 more seconds/cycles
            breakloop = true; % Break simulation loop after those 10 seconds/cycles
        end
        %}
    end
end % End while loop
if breakloop == false % Simulation did not reach eps_sync
    % Set t_eps to be the end of the simulation, t(end,:)
    t_eps = t(k,:);
    k_eps = k_r;% - sum(t > (t_eps-(10*(threshold/w_0))));
end
%c_avg = max(c_roll(k_eps:k));
c_avg = c_roll(k_eps);
% % Calculate average value of c(k) for last 10 seconds/cycle
% c_avg = 0;
% for id = k_eps:k-1
%     val = 0.5*(c(id)+c(id+1))*(t(id+1,:)-t(id,:));
%     c_avg = c_avg + val;
% end
% c_avg = c_avg / (t(k,:)-t(k_eps,:));
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
    case 'ERGrand' % Erdős–Rényi–Gilbert random network
        A = rand(N,N) + eye(N);
        if symmetric % birectional graph
            A = A - tril(A,-1); % remove lower triangular components
            A = A+A';
        end
        A = double(A<vector); % vector is probability of making an edge (0,1]
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