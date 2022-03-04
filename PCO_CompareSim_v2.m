% PCO Simulation Comparison - Version 2
% Based on PCO_Sim031020
% Compares performance of PCO synchronization PRFs/algorithms
% Focuses on behavior in non-ideal environments
% Last Modified: 11/16/2021

%% Simulation Parameters
cycles = 101; % Maximum number of evolution cycles
N = 20; % Number of PCOs (integer >= 2)
alpha = 0.99 * ones(N,1); % PCO coupling strength - (0,1]
%alpha = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6];
threshold = 2*pi; % Oscillator phase threshold
w_0 = 2*pi; % Nominal oscillator evolution frequency
D = 0.0; % Refractory Period
A_option = 'ERGrand';
%A = GetTopology(N,A_option);
%A = GetTopology(N,'custom',[2:5],true); % Custom topology
A = GetTopology(N,'ERGrand',0.15,true); % Random topology
out_degree = sum(A,2); % Outdegree of oscillators
out_degree_avg = sum(out_degree)/N; % Average network outdegree
in_degree = sum(A,1); % Indegree of oscillators
in_degree_avg = sum(out_degree)/N; % Average network indegree
eps_sync = 0.0001*threshold; % epsilon-synchronization threshold

%theta_start = theta; % Store initial oscillator states (for additional simulations)

area_1 = (0.5*threshold)^2;
c = [0.07 + 0.43*(exp(-0.23*((in_degree')-1))), 0.77 - 0.27*(exp(-0.13*((in_degree')-1)))]; % learned PRF parameters
area_2 = 0.5*threshold^2 * (c(:,1).^2 + (1-c(:,2)).^2 + (1-c(:,1)).*c(:,1) - c(:,1).*(1-c(:,2)).*(1-c(:,2))./(1-c(:,1)));

%global epsilon beta; % MS parameters
epsilon = 0.01*ones(N,1);
beta = 3*ones(N,1);
% Estimate value of epsilon that gives same "area" as GeneralPRF
%{
p_step = threshold/1000;
p = 0:p_step:threshold;
stateMS = log(1+(exp(beta)-1)*(p./threshold))./beta;
for i=1:10
    phaseMS = threshold.*(exp(beta.*(stateMS + epsilon))-1)./(exp(beta)-1);
    p_plus = min((phaseMS - p),(threshold - p));
    area_3 = p_step * (sum(p_plus,2) - 0.5*(p_plus(:,1)+p_plus(:,end)));
    epsilon = epsilon.*(area_2.*(1.*alpha)./area_3);
end
%}

%% Simulation
rng('default');
sample_size = 50; % Number of simulations to run/compare
w_var = 0.01; % oscillator frequency variation (percentage of w_0)
tau_min = 0.01; % minimum pulse delay (fraction of cycle)
tau_max = 0.04; % maximum pulse delay (fraction of cycle)
r_max = 0.5; % partition parameter for Javed2021
r_min = 0.5; % partition parameter for Javed2021
L = zeros(sample_size,1); % Initial containing arc of oscillator for each set of simulations
maxTime = zeros(sample_size,4); % Time to synchronize for each simulation
SyncAmount = zeros(sample_size,4); % Average containing arc value once network reaches below eps_sync.
sample = 0; % Sample number index
thetas = zeros(sample_size,N);
omegas = zeros(sample_size,N);
pulse_delays = zeros(N,N,sample_size);
r_parts = zeros(N,sample_size);
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
    rand_delay = tau_min*ones(N,N) + (tau_max-tau_min)*rand(N,N); % Random matrix
    pulse_delays(:,:,i) = 0.5 * (rand_delay + rand_delay'); % (symmetric) matrix of pulse propogation delays
    r_parts(:,i) = r_min*ones(N,1) + (r_max - r_min)*rand(N,1);
    sample = sample + 1; % Increment sample counter
    L(sample) = Lambda_0;
end
start_time = tic; % For timing the simulations
for sample = 1:sample_size
    theta = thetas(sample,:);
    omega = omegas(sample,:);
    pulse_delay = pulse_delays(:,:,sample);
    r_part = r_parts(:,sample);
    %
    [time1, state1, C_Arc1, t_sync1, C_avg1, C_roll1] = PCO_Sync_Sim(@StandardPRF, theta, w_0, omega, pulse_delay, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    [time2, state2, C_Arc2, t_sync2, C_avg2, C_roll2] = PCO_Sync_Sim(@GeneralPRF, theta, w_0, omega, pulse_delay, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    %
    [time3, state3, C_Arc3, t_sync3, C_avg3, C_roll3] = PCO_Sync_Sim(@KlinglmayrPRF, theta, w_0, omega, pulse_delay, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    %
    [time4, state4, C_Arc4, t_sync4, C_avg4, C_roll4] = PCO_Sync_Sim(@MSPRF, theta, w_0, omega, pulse_delay, cycles, A, threshold, D, 1.0, 1.0, ...
                                            alpha, epsilon, beta, tau_min, tau_max, c(:,1), c(:,2), r_part, eps_sync);
    %}
    maxTime(sample,:) = [t_sync1 t_sync2 t_sync3 t_sync4];
    SyncAmount(sample,:) = [C_avg1 C_avg2 C_avg3 C_avg4];
    % Calculate percentage complete and display
    DisplayPercentageComplete(sample, sample_size, toc(start_time), 2);
end

%%
%sync_bool = not(maxTime < cycles);
sync_bool = SyncAmount < 0.5*threshold;
[~,min_index] = min(maxTime,[],2);
%disp(maxTime)
disp_data = zeros(5,4); total_sync = zeros(1,4); fast_sync = zeros(1,4);
for i = 1:4
disp_data(1,i) = median(maxTime(sync_bool(:,i),i)); %disp([median(maxTime(sync_bool(:,1),1)), median(maxTime(sync_bool(:,2),2)), median(maxTime(sync_bool(:,3),3))])
disp_data(2,i) = mean(maxTime(sync_bool(:,i),i)); %disp([mean(maxTime(sync_bool(:,1),1)), mean(maxTime(sync_bool(:,2),2)), mean(maxTime(sync_bool(:,3),3))])
disp_data(3,i) = std(maxTime(sync_bool(:,i),i)); %disp([std(maxTime(sync_bool(:,1),1)), std(maxTime(sync_bool(:,2),2)), std(maxTime(sync_bool(:,3),3))])
disp_data(4,i) = mean(SyncAmount(sync_bool(:,i),i)*100/threshold); % Average percent of cycle length
disp_data(5,i) = std(SyncAmount(sync_bool(:,i),i)*100/threshold); % std. dev of cycle length in percent
total_sync(i) = sum(sync_bool(:,i)); %disp(sum(sync_bool))
fast_sync(i) = sum(min_index == i); %disp([sum(min_index == 1), sum(min_index == 2), sum(min_index == 3)])
end
disp(disp_data)
disp(total_sync)
disp(fast_sync)
%disp([mean(SyncAmount(:,1)), mean(SyncAmount(:,2)), mean(SyncAmount(:,3)), mean(SyncAmount(:,4))])

%% Plots
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
plot(time1,C_Arc1, 'k-', 'Linewidth', 1.5) % StandardPRF
plot(time2,C_Arc2, 'r-', 'Linewidth', 1.5) % GeneralPRF
plot(time3,C_Arc3, 'b-', 'Linewidth', 1.5) % Klinglmayr2017
plot(time4,C_Arc4, '-', 'Linewidth', 1.5) % Javed2021
hold off
grid on
xlabelh = xlabel('Time (s)');
ylabelh = ylabel('Containing Arc (\Lambda)');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Klinglmayr (2017)', 'Javed (2021)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northeast')
axis([0,max(maxTime(end,:)),0,threshold*(N-1)/N])
set(gca, 'yscale','log')
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
plot(L,maxTime(:,4), '*')
hold off
grid on
xlabelh = xlabel('Initial Containing Arc (\Lambda)');
ylabelh = ylabel('Cycles to \epsilon-synchronization');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF', 'Klinglmayr (2017)', 'Javed (2021)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northwest')
axis([min(L),threshold*(N-1)/N,0,cycles-1])
set(gca, 'yscale','log')
%}

%% Hardcoded Data Collection
alpha_values = [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02, 0.01, 0.005];
sync_prob_Opt = [1, 0.9985, 0.9942, 0.9887, 0.9813, 0.9685, 0.9495, 0.9363, 0.9288, 0.9291, 0.9287, 0.9484, 0.9598, 0.9619];
sync_prob_RL = [1, 1, 1, 1, 0.9869, 0.9755, 0.9685, 0.9635, 0.9645, 0.9675, 0.9698, 0.9730, 0.9761, 0.9766];

% Figure 4: Compare Prob. Sync
%{
fig4 = figure(4);
clf
set(fig4, 'Position', [55, 50, 500, 440])
hold on
plot(alpha_values, sync_prob_Opt, 'kx-', 'LineWidth', 1.5)
plot(alpha_values, sync_prob_RL, 'r*-', 'LineWidth', 1.5)
hold off
grid on
xlabelh = xlabel('Coupling Strength, l');
ylabelh = ylabel('Probability of Synchronization');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF', 'Learned PRF');%, 'Klinglmayr (2017)', 'Javed (2021)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northwest')
%axis([min(L),threshold*(N-1)/N,0,cycles-1])
%set(gca, 'yscale','log')
%}

% Quartiles of data
maxTime_q = zeros(4,3);
SyncAmount_q = zeros(4,3);
for i = 1:4
    maxTime_q(i,:) = quantile(maxTime(sync_bool(:,i),i),[0.25, 0.5, 0.75]);
    SyncAmount_q(i,:) = quantile(SyncAmount(sync_bool(:,i),i),[0.25, 0.5, 0.75]);
    %mid_y(i) = disp_data(4,i); % Average
    %mid_x(i) = disp_data(2,i);
    mid_y(i) = SyncAmount_q(i,2); % Median
    mid_x(i) = maxTime_q(i,2);
    err_yh(i) = SyncAmount_q(i,3)-mid_y(i);
    err_yl(i) = mid_y(i)-SyncAmount_q(i,1);
    err_xh(i) = maxTime_q(i,3)-mid_x(i);
    err_xl(i) = mid_x(i)-maxTime_q(i,1);
end
markers = {'kx', 'rx', 'bx', 'x'};
% Figure 5: Comparison
%{
fig5 = figure(5);
clf
set(fig5, 'Position', [555, 50, 500, 400])
hold on
for i = 1:4
    errorbar(mid_x(i), mid_y(i), err_yl(i),err_yh(i), err_xl(i),err_xh(i), markers{i}, 'Linewidth', 1)
end
hold off
grid on
xlabelh = xlabel('Average cycles to steady-state Containing Arc');
ylabelh = ylabel('Average value of steady-state Containing Arc');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF (l=0.5)', 'Learned PRF (l=0.5)', 'Klinglmayr (2017)', 'Nishimura (2011)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northeast')
titleh = title('N = 6 Oscillators, All-to-All topology');
%axis([0,25, 0,0.09])
%set(gca, 'yscale','log')
%}

% Figure 6: Rolling average of Containing Arc
%{
fig6 = figure(6);
clf
set(fig6, 'Position', [55, 400, 500, 240])
hold on
plot(time1(1:length(C_roll1)), C_roll1, 'kx-', 'LineWidth', 1)
plot(time2(1:length(C_roll2)), C_roll2, 'ro-', 'LineWidth', 1)
plot(time3(1:length(C_roll3)), C_roll3, 'b+-', 'LineWidth', 1)
plot(time4(1:length(C_roll4)), C_roll4, '*-', 'LineWidth', 1)
hold off
grid on
xlabelh = xlabel('Time (cycles)');
ylabelh = ylabel('Rolling Average of Containing Arc');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF (l=0.5)', 'Learned PRF (l=0.5)', 'Klinglmayr (2017)', 'Javed (2021)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'northeast')
%titleh = title('N = 6 Oscillators, All-to-All topology');
%axis([0,25, 0,0.09])
%set(gca, 'yscale','log')
%}

%% PRF plots
p_steps = 1000;
p = linspace(0,threshold,p_steps+1);
StandF = zeros(size(p));
RLF = zeros(size(p));
KlingF = zeros(size(p));
%JavedF = zeros(N,length(p));
NishiF = zeros(size(p));
r_part = r_min*ones(N,1) + (r_max - r_min)*rand(N,1);
params = [alpha, epsilon, beta, tau_min*ones(N,1), tau_max*ones(N,1), c(:,1), c(:,2), r_part];
for i = 1:length(p)
    StandF(i) = StandardPRF(p(i),1,threshold,D,params);
    RLF(i) = GeneralPRF(p(i),1,threshold,D,params);
    KlingF(i) = KlinglmayrPRF(p(i),1,threshold,D,params);
    %for j = 1:N
    %    JavedF(j,i) = JavedPRF(p(i),j,threshold,D,params);
    %end
    NishiF(i) = NishimuraPRF(p(i),1,threshold,D,params);
end
KlingF_mod = mod(KlingF + (threshold/2), threshold) - (threshold/2);
KlingF_area = (threshold/p_steps) * (sum(abs(KlingF_mod),2) - 0.5*(abs(KlingF_mod(:,1))+abs(KlingF_mod(:,end))));
%disp(KlingF_area/area_1)
NishiF_area = (threshold/p_steps) * (sum(abs(NishiF),2) - 0.5*(abs(NishiF(:,1))+abs(NishiF(:,end))));
%disp(NishiF_area/area_1)

% Figure 9: PRF comparison
%
fig9 = figure(9);
clf
set(fig9, 'Position', [55, 50, 350, 440])
hold on
plot(p, StandF, 'k-', 'LineWidth', 1)
plot(p, RLF, 'r--', 'LineWidth', 1)
plot(p, KlingF, 'b-', 'LineWidth', 1)
%plot(p, KlingF_mod, 'c-', 'LineWidth', 1)
plot(p, NishiF, '-.', 'LineWidth', 1)
hold off
grid on
xlabelh = xlabel('Oscillator Phase');
ylabelh = ylabel('Phase Response Function');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
legendh = legend('Standard PRF (l=0.5)', 'Learned PRF (l=0.5)', 'Klinglmayr (2017)', 'Nishimura (2011)');
set(legendh,'Fontname','Times New Roman', 'Fontsize',12, 'Location', 'southwest')
axis([0,threshold, -threshold,threshold/2])
xticks([0 (pi/2) pi (3*pi/2) 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
yticks([-2*pi -3*pi/2 -pi (-pi/2) 0 (pi/2) pi 3*pi/2 2*pi])
yticklabels({'-2\pi', '-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2', '2\pi'})

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
    y = (phase_plus - phase);
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
    y = H - phase;
end
end

% Javed 2021 PRF for PCO Synchronization
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
%   params = vector; [alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r]
%       Note: not all parameters are used in each function.
% Returns:
%   y = float; the value of the PRF at phase
function [y] = JavedPRF(phase, index, threshold, refractory, params)
r = params(index, 8);
if phase > r*threshold
    y = threshold - phase;
elseif phase > refractory
    y = -1 * phase;
else
    y = 0.0;
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
kappa = 0.01*threshold;
if phase > r*threshold % Something excitory (positive)
    y = (threshold/pi)*sin((pi/threshold)*phase);
    %y = params(index,1)*(threshold - phase);
elseif phase > tau + kappa % Something inhibitory (negative) enough
    y = min(-(tau+kappa), -(threshold/pi)*sin((pi/threshold)*phase));
    %y = min(-(tau+kappa), -params(index,1)*phase);
elseif phase > refractory
    y = -1 * phase;
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

function [t, s, c, t_eps, c_avg, c_roll] = PCO_Sync_Sim(PRF, theta, w_0, omega, pulse_delay, cycles, A, threshold, D, p_send, p_receive, alpha, epsilon, beta, tau_min, tau_max, RL_c1, RL_c2, r, eps_sync)
N = length(theta); % Number of oscillators in the network
PRF_params = [alpha, epsilon, beta, tau_min*ones(N,1), tau_max*ones(N,1), RL_c1, RL_c2, r];
p_s = p_send; % probability of sending a pulse when resetting (Klinglmayr)
p_r = p_receive; % probability of receiving a sent pulse (Javel)

roll_window = 5;
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
                        %delay{i} = [delay{i}, pulse_delay(index_f,i)]; % Fixed delay
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
        if c(k_r) < (threshold/10) && max(c_roll(k_start:k_r)) - min(c_roll(k_start:k_r)) < eps_sync && breakloop == false
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