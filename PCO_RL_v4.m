% Reinforcement Learning Algorithm for PCO Network
% By: Timothy Anglea
% Last Modified: 7/20/2020
% Determination of Reward value after each firing using estimated state
% Parameterized Sarsa batch update
% Heterogeneous state-action values - Each oscillator learns its own optimal PRF

%% Reinforcement Learning Parameters
P = 100; % PRF parameter size (minus 1)
N = 6; % Network size (number of oscillators)
threshold = 2*pi; % Oscillator phase threshold

Action_factor = 2; % Interger multiple
% Action_space stepsize = 1/(Action_factor*P);
num_Action = Action_factor*(P) + 1; % number of action steps for each state
%rng('default');
Q_option = 2;
if exist('Q','var') == 0 || exist('Q_start','var') == 0
    Q_option = 2; Q = 'none'; Q_start = 'none';
end
[Q,P,num_Action] = GetQ(Q_option,N,P,num_Action,Q,Q_start);
Q_start = Q;
State_positions = linspace(0,threshold,P+1)'; % Parameterization of state positions
Action_space = GetActionSpace(3,State_positions,num_Action,threshold);
% Action_space(s,:) are possible actions (phase changes) for an oscillator with state "s"

%% PCO Simulation (modified from PCO_Sim031020.m)
%% Simulation Parameters
cycles = 15; % Number of evolution cycles
l = 1.0; % PCO coupling strength - (0,1]
%D = 0.0; % Refractory Period
A_option = 'all';
%a_block = GetTopology(N/3,A_option);
%A = blkdiag(a_block,a_block,a_block);%,a_block,a_block);
A = GetTopology(N,A_option);
%A = GetTopology(N,'custom',2:5,false); % Custom topology
%A = GetTopology(N,'ERGrand',0.15,true); % Random topology
%A = GetTopology(N,'alldel',0,true); % Every indegree topology
out_degree = sum(A,2); % Outdegree of oscillators
out_degree_avg = sum(out_degree)/N; % Average network outdegree
in_degree = sum(A,1) % Indegree of oscillators
in_degree_avg = sum(in_degree)/N; % Average network indegree
w_0 = 2*pi; % Nominal oscillator evolution frequency
% Nominal (simulation) time to complete simulation (in seconds): cycles * (threshold/w_0);
T = 5000; % Number of episodes
%alpha_SAR = 1/1000; %0.01/in_degree_avg; % SARSA learning rate
alpha_SAR_list = [0.01, 0.01, 0.005, 0.005, 0.002, 0.002, 0.001, 0.001, 0.0005, 0.0005];
gamma_SAR = 0.5; % SARSA discount rate
%beta = 8; % Policy selection factor
beta_list = [2,2,3,3,3,3,4,4,4,5];
cycle_fraction = 1.0; % Fraction of cycle for initial phases (0,1]
freq_var = 0.05; % Amount of oscillator frequency variation (non-negative)
delay_const = 0.0; % Amount of constant delay for each oscillator (non-negative)
delay_var = 0.1; % Amount of delay variation for each oscillator (non-negative)
%%
Rewards = zeros(T,N);
Q_time = zeros([size(Q),100+1]); % Q-values after each episode
Q_time(:,:,:,1) = Q; % Initial Q-values for first episode
Q_step = round(T/100);
Qt_vec = 1:Q_step:T+1;
start_time = tic; % For timing the simulation
for t = 1:T
    ti = ceil(t/(T/10));
    alpha_SAR = alpha_SAR_list(ti);
    beta = beta_list(ti);
    %{
    if t <= 0.4*T
        %alpha_SAR = 0.01;
        beta = 2; % Probability Choice
        %beta = 0.9; % Epsilon-greedy
    elseif t <= 0.7*T
        %alpha_SAR = 0.01;
        beta = 3; % Probability Choice
        %beta = 0.5; % Epsilon-greedy
    else % t > 0.7*T
        %alpha_SAR = 0.01;
        beta = 5; % Probability Choice
        %beta = 0.3; % Epsilon-greedy
    end
    %}
    %%Start of PCO simulation
    theta = cycle_fraction*rand(1,N)*threshold; % vector of initial node phases [0,1]
    omega = w_0 * (ones(1,N) + freq_var*(2*(rand(1,N)-0.5))); % vector of oscillator frequencies
    %rand_delay = (delay_const)*ones(N,N) + (delay_var)*rand(N,N); % Random matrix
    %pulse_delay = (0.5 * (rand_delay + rand_delay')); % (symmetric) matrix of pulse propogation delays
    %pulse_delay = rand_delay - tril(rand_delay,-1) + triu(rand_delay,1)'; % Actually uniformly random
    % Find PRC values for use in episode
    PRC = GetPRC(Q, Action_space, beta);
    
    % Simulate the network to determine experience for the episode
    T_episode = cycles * (threshold/w_0);
    [~, ~, ~, SA_pairs, R] = PCO_Sync_RL(T_episode, omega, threshold, theta, A, l, delay_const, delay_var, State_positions, PRC);
    
    %%After end of episode, use SA_pairs and R to update Q (SARSA, etc.)
    Rewards(t,:) = sum(R); % Total reward gained over episode for each oscillator.
    % Can be useful measure to observe learning over time.
    
    % Update state-action values based on previous episode
    Q = QUpdateSARSA(Q, SA_pairs, R, alpha_SAR, gamma_SAR, State_positions, Action_space, PRC, threshold);
    % Record Q-values periodically after a training episode
    for i = 1:length(Qt_vec)
        if t+1 < Qt_vec(i)
            break
        end
        if t+1 == Qt_vec(i)
            Q_time(:,:,:,i) = Q; % Record updated Q-values
            break
        end
    end
    % Calculate percentage complete and display
    DisplayPercentageComplete(t,T,toc(start_time),1);
end

%% Plots and Figures

% Figure 1: Optimal PRF based on Q-values
% Find PRF values for each State position parameter
PRF1 = zeros(size(PRC));
for n = 1:N
    for i = 1:P+1
        [Qin_sort, A_opt_indices] = sort(Q(i,:,n),'descend');
        PRF1(i,n) = Action_space(i,A_opt_indices(1));
    end
end

PRF_avg = mean(PRF1,2);
PRF_min = min(PRF1,[],2);
PRF_max = max(PRF1,[],2);
% Hard-coded learned PRF parameters
c = [0.095 + 0.405*(exp(-0.29*((in_degree')-1))), 0.755 - 0.255*(exp(-0.14*((in_degree')-1)))];
c1 = c(1,1)*threshold; % just use first set of parameters
c2 = c(1,2)*threshold;
cx = [0,c1,c2,c2,threshold];
cy = [0,-c1,(c1/(threshold-c1))*(c2-threshold),threshold-c2,0];
black = [0,0,0]; grey = [0.55,0.55,0.55]; magenta = [0.9,0.25,0.9];
%
fig1h = figure(1);
clf
set(fig1h, 'Position', [100, 290, 500, 350])
hold on
%plot(State_positions, PRF1, '-', 'Color', [0,0.7,0], 'Linewidth', 0.5)
plot(State_positions, PRF_avg, 'x-', 'Color', black, 'Linewidth', 1.5)
plot(cx(1:3),cy(1:3),'-', 'Color', magenta, 'Linewidth', 1.5)
plot(cx(3:4),cy(3:4),'--', 'Color', magenta, 'Linewidth', 1.5)
plot(cx(4:5),cy(4:5),'-', 'Color', magenta, 'Linewidth', 1.5)
plot(State_positions, PRF_max, ':', 'Color', grey, 'Linewidth', 1.5)
plot(State_positions, PRF_min, ':', 'Color', grey, 'Linewidth', 1.5)
hold off
grid on
xlabelh = xlabel('Phase, $\phi$', 'interpreter', 'latex');
ylabelh = ylabel('Phase Update, $F(\phi)$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,threshold,-0.53*threshold,0.53*threshold])
xticks([0 (pi/2) pi (3*pi/2) 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
%yticks([-pi (-pi/2) 0 (pi/2) pi])
%yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
yticks([-3*pi/2 -pi (-pi/2) 0 (pi/2) pi 3*pi/2])
yticklabels({'-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2'})

legendh = legend('Avg. Learned Policy', 'Proposed Model Fit');
set(legendh,'Fontname','Times New Roman', 'Fontsize',11, 'Location', 'northwest')
%}

o = 1; % Oscillator number (Used for Figure 2 & 3)
% Figure 2: Q-values !!! Not Testeded with multiple Q-value sets !!!
%{
Q_sort = sort(reshape(Q(:,:,o),[(P+1)*num_Action,1]));
offset = round((P+1)*num_Action*0.001);
max_Q = Q_sort(end-offset);
min_Q = Q_sort(1+offset);
%max_Q = max(max(Q));
%min_Q = min(min(Q));
fig2h = figure(2);
clf
set(fig2h, 'Position', [310, 200, 500, 400])
imagesc(Q(:,:,o))
colorbar
colormap parula
caxis([min_Q,max_Q]);
xlabelh = xlabel('Action Index');
ylabelh = ylabel('State Index');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
%}

% Figure 3: Q-values plotted over Action space
%{
%max_Q = max(max(Q(2:end,:)));
%min_Q = min(min(Q(2:end,:)));
Q_vec = reshape(Q(:,:,o),[],1);
A_vec = reshape(Action_space,[],1);
S_vec = repmat(State_positions,num_Action,1);
fig3h = figure(3);
clf
set(fig3h, 'Position', [900, 230, 360, 400])
hold on
scatter(S_vec,A_vec,7,Q_vec,'s','filled')
%{
for s = 1:P
    x_val = State_positions([s,s,s+1,s+1,s]);
    for a = 1:num_Action
        if a == 1 % Color lower triangles
            for i = 1:Action_factor
                y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a,a+i-1,a+i,a]));
                colors = diag(Q([s,s,s+1,s+1,s],[a,a,a+i-1,a+i,a],o));
                fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
            end
        elseif a > num_Action - Action_factor % Color upper triangles
            y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a-1,num_Action,num_Action,a]));
            colors = diag(Q([s,s,s+1,s+1,s],[a,a-1,num_Action,num_Action,a],o));
            fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
        else % Color rectangle
            y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a-1,a-1+Action_factor,a+Action_factor,a]));
            colors = diag(Q([s,s,s+1,s+1,s],[a,a-1,a-1+Action_factor,a+Action_factor,a],o));
            fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
        end
    end
    %drawnow
end
%}
%plot(State_positions, PRC, 'o-', 'Linewidth', 1.5, 'Color', [1,0,0])
for n = 1:1 %1:N
plot(State_positions, PRF1(:,o), 'x-', 'Linewidth', 1.5, 'Color', (n-1)*(ones(1,3)/N))
end
grid on
hold off
colorbar
caxis([min_Q,max_Q]);
xlabelh = xlabel('Phase, $\phi$', 'interpreter', 'latex');
ylabelh = ylabel('Phase Update, $F(\phi)$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,threshold,-threshold,threshold])
xticks([0 (pi/2) pi (3*pi/2) 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
yticks([-2*pi -3*pi/2 -pi (-pi/2) 0 (pi/2) pi 3*pi/2 2*pi])
yticklabels({'-2\pi', '-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
%}

% Figure 4: Total Reward earned in an episode over iterations
%{
fig4h = figure(4);
clf
%set(fig4h, 'Position', [1310, 300, 500, 300])
plot(1:T,Rewards, '.', 'Linewidth', 1.0)
grid on
xlabelh = xlabel('Episode');
ylabelh = ylabel('Total Reward');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,T,-6.3,6.3])
%}

% Figure 5: Evolution of containing arc in an episode
%{
fig5h = figure(5);
clf
set(fig5h, 'Position', [110, 300, 500, 300])
hold on
for i = 1:10
    % Simulate the network to evaluate performance of optimal PRF
    theta = cycle_fraction*rand(N,1)*threshold; % vector of initial node phases [0,1]
    omega = w_0 * ones(N,1) + (freq_var*threshold)*(rand(N,1)-0.5); % vector of oscillator frequencies
    %rand_delay = (delay_const)*ones(N,N) + (delay_var)*rand(N,N); % Random matrix
    %pulse_delay = (0.5 * (rand_delay + rand_delay')); % (symmetric) matrix of pulse propogation delays
    T_episode = cycles * (threshold/w_0);
    [time, state, c, ~, ~] = PCO_Sync_RL(T_episode, omega, threshold, theta, A, l, delay_const, delay_var, State_positions, PRF1);
    plot(time, c, '-', 'Linewidth', 1.5)
end
hold off
grid on
axis([0,cycles,0,threshold*(N-1)/N])
xlabelh = xlabel('Time');
ylabelh = ylabel('Containing Arc');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
%}

% Figure 6: Evolution of oscillator phases in an episode
%{
fig6h = figure(6);
clf
%set(fig5h, 'Position', [1310, 300, 500, 300])
plot(time, state, 'x-', 'Linewidth', 1.5)
grid on
xlabelh = xlabel('Time');
ylabelh = ylabel('Oscillator Phases');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
%}

%{
o = 1; % Oscillator number for plot
fig7h = figure(7);
set(fig7h, 'Position', [310, 150, 360, 400])
for t = 1:100+1
    clf
    Q_time_sort = sort(reshape(Q_time(:,:,o,t),[(P+1)*num_Action,1]));
    offset = round((P+1)*num_Action*0.002);
    max_Q_time = Q_time_sort(end-offset);
    min_Q_time = Q_time_sort(1+offset);
    PRF_opt = zeros(P+1,1);
    for p = 1:P+1
        [Qin_sort, A_opt_indices] = sort(Q_time(p,:,o,t),'descend');
        PRF_opt(p) = Action_space(p,A_opt_indices(1));
    end
    
    Q_vec = reshape(Q_time(:,:,o,t),[],1);
    A_vec = reshape(Action_space,[],1);
    S_vec = repmat(State_positions,num_Action,1);
    hold on
    scatter(S_vec,A_vec,7,Q_vec,'s','filled')
    %{
    for s = 1:P
        x_val = State_positions([s,s,s+1,s+1,s]);
        for a = 1:num_Action
            if a == 1 % Color lower triangles
                for i = 1:Action_factor
                    y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a,a+i-1,a+i,a]));
                    colors = diag(Q_time([s,s,s+1,s+1,s],[a,a,a+i-1,a+i,a],1,t));
                    fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
                end
            elseif a > num_Action - Action_factor % Color upper triangles
                y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a-1,num_Action,num_Action,a]));
                colors = diag(Q_time([s,s,s+1,s+1,s],[a,a-1,num_Action,num_Action,a],1,t));
                fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
            else % Color rectangle
                y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a-1,a-1+Action_factor,a+Action_factor,a]));
                colors = diag(Q_time([s,s,s+1,s+1,s],[a,a-1,a-1+Action_factor,a+Action_factor,a],1,t));
                fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
            end
        end
    end
    %}
    plot(State_positions, PRF_opt, 'x-', 'Linewidth', 1.5, 'Color', [0,0,0])
    hold off
    grid on
    colorbar
    caxis([min_Q_time,max_Q_time]);
    titleh = title(['Q-Values after Episode ',num2str(Qt_vec(t)-1)]);
    set(titleh,'Fontname','Times New Roman', 'Fontsize',14)
    xlabelh = xlabel('Phase, \phi');
    ylabelh = ylabel('Phase Update, F(\phi)');
    set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
    set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
    axis([0,threshold,-threshold,threshold])
    xticks([0 (pi/2) pi (3*pi/2) 2*pi])
    xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
    yticks([-2*pi -3*pi/2 -pi (-pi/2) 0 (pi/2) pi 3*pi/2 2*pi])
    yticklabels({'-2\pi', '-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
    drawnow
end
%}

% Test for RL schematic
%
fig8h = figure(8);
clf
set(fig8h, 'Position', [650, 290, 500, 250])
hold on
%{
o=1;
for s = 1:P
    x_val = State_positions([s,s,s+1,s+1,s]);
    for a = 1:num_Action
        if a == 1 % Color lower triangles
            for i = 1:Action_factor
                y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a,a+i-1,a+i,a]));
                colors = diag(Q([s,s,s+1,s+1,s],[a,a,a+i-1,a+i,a],o));
                fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
            end
        elseif a > num_Action - Action_factor % Color upper triangles
            y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a-1,num_Action,num_Action,a]));
            colors = diag(Q([s,s,s+1,s+1,s],[a,a-1,num_Action,num_Action,a],o));
            fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
        else % Color rectangle
            y_val = diag(Action_space([s,s,s+1,s+1,s],[a,a-1,a-1+Action_factor,a+Action_factor,a]));
            colors = diag(Q([s,s,s+1,s+1,s],[a,a-1,a-1+Action_factor,a+Action_factor,a],o));
            fill(x_val,y_val,colors, 'LineStyle', 'none') % Plot
        end
    end
    %drawnow
end
%}
plot(State_positions, PRF1(:,1), 'x:', 'Color', [0,0,0], 'Linewidth', 2, 'MarkerSize', 14)
hold off
%grid on
%xlabelh = xlabel('Phase, $\phi$', 'interpreter', 'latex');
%ylabelh = ylabel('Phase Update, $F(\phi)$', 'interpreter', 'latex');
%set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
%set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,threshold,-0.25*threshold,0.32*threshold])
set(gca, 'box', 'off', 'xcolor', 'none', 'ycolor', 'none', 'Position', [0.02, 0.01, 0.96,0.98])

%xticks([0 (pi/2) pi (3*pi/2) 2*pi])
%xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
%yticks([-pi (-pi/2) 0 (pi/2) pi])
%yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
%yticks([-3*pi/2 -pi (-pi/2) 0 (pi/2) pi 3*pi/2])
%yticklabels({'-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2'})
%}

%% Functions
% Create state-action value matrix, Q
function [Q,P,num_Action] = GetQ(option,N,P,num_Action,Q,Q_start)
switch option
    case 1 % Random
        Q = rand(P+1,num_Action,N)-0.5; % Randomized initial value of state-action pairs
    case 2 % Equal
        Q = (0.0)*ones(P+1,num_Action,N); % Constant initial value of state-action pairs
    case 3 % Continue with last Q
        disp("Using previous final values of Q")
        [P_prime, num_Action, ~] = size(Q);
        P = P_prime - 1;
    case 4 % Double size of Q
        disp("Doubling size of Q")
        Q = ScaleQ(Q); % !!!Not sure this works yet!!!
        [P_prime, num_Action, ~] = size(Q);
        P = P_prime - 1;
    otherwise % Use last initial state-action value pairs
        disp("Using previous initial values of Q")
        Q = Q_start;
        [P_prime, num_Action, ~] = size(Q);
        P = P_prime - 1;
end
end

% Generate action values for each state
function [Acts] = GetActionSpace(option,State_positions,num_Action,threshold)
P_num = length(State_positions);
Acts = zeros(P_num,num_Action);
% Acts(s,:) are possible actions (phase changes) for an oscillator with state "s"
switch option
    case 1 % Restrict actions to delay-advance
        for i = 1:P_num
            if State_positions(i) < 0.5*threshold
                min_action = 0 - State_positions(i);
                max_action = 0;
            elseif State_positions(i) == 0.5*threshold
                min_action = 0 - State_positions(i);
                max_action = threshold - State_positions(i);
            else % State_positions(i) > 0.5*threshold
                min_action = 0;
                max_action = threshold - State_positions(i);
            end
            Acts(i,:) = linspace(min_action,max_action,num_Action);
        end
    case 2 % Cap actions to half-cycle adjustments
        for i = 1:P_num
            if State_positions(i) < 0.5*threshold
                min_action = 0 - State_positions(i);
                max_action = 0.5*threshold;
            elseif State_positions(i) == 0.5*threshold
                min_action = 0 - State_positions(i);
                max_action = threshold - State_positions(i);
            else % State_positions(i) > 0.5*threshold
                min_action = -0.5*threshold;
                max_action = threshold - State_positions(i);
            end
            Acts(i,:) = linspace(min_action,max_action,num_Action);
        end
    otherwise % Allow all possible actions
        for i = 1:P_num
            min_action = 0 - State_positions(i);
            max_action = threshold - State_positions(i);
            Acts(i,:) = linspace(min_action,max_action,num_Action);
        end
end
end

% Generate PRC for an episode, based on current state-action values
function [PRC] = GetPRC(Q, Action_space, beta)
PRC = zeros(size(squeeze(Q(:,1,:))));
for n = 1:length(PRC(1,:)) % For each oscillator
    for i = 1:length(PRC(:,1)) % For each state
        % Pick PRC (action) based on current optimal (epsilon-greedy) !!!Not tested!!!
        %{
        if rand(1) > beta
            A_opt_index_possible = find(Q(i,:,n) == max(Q(i,:,n)));
            A_opt_index = A_opt_index_possible(randsample(length(A_opt_index_possible),1));
            PRC(i,n) = Action_space(i,A_opt_index);
        else
            PRC(i,n) = Action_space(i,randi([1,length(Action_space)],1));
        end
        %}
        % Pick PRC (action) based on exponential function probability
        %
        P_value = exp(beta * (Q(i,:,n) - mean(Q(i,:,n))));
        P_pick = P_value./sum(P_value); % Probability of picking an action for the PRC
        rand_num = rand(1);
        for j = 1:length(P_pick)
            rand_num = rand_num - P_pick(j);
            if rand_num < 0
                A_index = j;
                break
            end
        end
        PRC(i,n) = Action_space(i,A_index);
        %}
    end
end
end

% Retreive action based on PRC
function [Action] = GetAction(theta, PRC, State_positions)
% Find current state based on state parameters
% State parameter below theta
[state_low, state_low_index] = max(State_positions(State_positions < theta));
% Catch condition if theta == 0
if isempty(state_low) % Indicates that theta == 0
    state_low_index = 1;
    state_low = State_positions(state_low_index);
end
% State parameter above theta(index_d)
state_high_index = state_low_index + 1;
state_high = State_positions(state_high_index);
% Choose phase change (i.e., action) based on policy;
% "Closeness" of theta to low state parameter
prob = (state_high - theta)/(state_high - state_low);
% Could choose PRC for episode in advance, then
A_low = PRC(state_low_index);
A_high = PRC(state_high_index);
Action = (prob*A_low) + ((1-prob)*A_high); % Combo action
end

% Compute containing arc, Lambda, of the PCO network
% Parameters:
%   theta = float; vector of oscillator phases
%   thresh = float; oscillator threshold
% Returns:
%   Lambda = float; containing arc of theta
function [Lambda] = ContainingArc(theta, thresh)
constant = ones(size(theta))*theta(1);
if isequal(theta, constant) % If the phases are all equal
    Lambda = 0; % The containing arc is zero
else
    theta_sorted = sort(theta); % sort the phases
    theta_shift = circshift(theta_sorted,1); % Shift phases over by 1
    v = mod(theta_sorted - theta_shift, thresh); % Phase difference between consecutive phases
    Lambda = thresh - max(v); % The containing arc is 1 minus the max difference
end
end

function [time, state, c, SA_pairs, R] = PCO_Sync_RL(T_eps, omega, threshold,...
    theta, A, l, d_c, d_v, State_positions, PRC)
% Simulation Parameters
%D = 0; % Refractory Period (to be incorporated)
%P_drop = 0; % Probability of pulse drop (to be incorporated)
N = length(theta); % Number of oscillators in the network
%in_degree = sum(A,1); % In-degree of oscillators in the network (unused)
%out_degree = sum(A,2); % Out-degree of oscillators in the network (unused)

% Start of PCO simulation
k = 1; % Simulation Iteration Index
k_SAR = ones(size(theta)); % Number of times oscillator i has fired a pulse
%m_SAR = ones(size(theta)); % Number of times (x2) oscillator i has received a pulse
% Store initial state
time = zeros(1,1); % Vector of times when data is stored
state = zeros(1,N); % Reset state values
state(k,:) = theta; % State of network at time(k)
c = zeros(1,1); % Reset arc values
c(k) = ContainingArc(state(k,:), threshold); % Containing Arc of network at time(k)
delta_matrix = cell(1,N); % Current estimate of phase differences from perspective of each oscillator
    % delta_theta_est{i} is oscillator i's estimate of the phase
    % differences to known oscillators from the previous cycle
m = zeros(1,N); % Number of oscillators heard since its own firing
c_est = zeros(1,N); % Containing arc of estimated network state for oscillators
    % c_est(k,i) is the estimate of the containing arc from oscillator i's perspective at timestep k
%{
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
    c_est(k,osc_order(i)) = ContainingArc([delta_matrix{osc_order(i)}, 0],threshold);
end
%}
%r_Lambda = c(1)*ones(size(theta)); % Actual value of containing arc; Used to compute reward for each oscillator
r_Lambda = c_est(k,:); % Estimated values of containing arc; Used to compute reward for each oscillator
delay = cell(size(theta)); % Empty cell array for delay timer lists
% SA_pairs = zeros(N,2,K); N = oscillators; 2 = state-action pair; K = cycle number;
SA_pairs = zeros([N,2,1]);
%Action = 0; % Clear action;
R = zeros(1,N); % Reward value computed for state-action pairs

% Main Loop: PCO Episode
while time(k) < T_eps % cycles * (threshold/w_0)
    %%Determine time to next event (responed to received pulse or reach threshold)
    % First, find smallest time to next received pulse
    receivetime = threshold/(min(omega)) + eps; % Largest possible time between firings
    index_d = 0; % index of oscillator with the smallest delay
    j = 0; % index of smallest delay for oscillator index_d
    for i = 1:1:length(delay) % Check each oscillator's delays
        if isempty(delay{i}) == false % If there is a delay for the oscillator
            if min(delay{i}) < receivetime % Is the delay the smallest so far?
                % If so, store that minimum time; the index of that timer
                % for the oscillator, and the index of the oscillator
                [receivetime, j] = min(delay{i});
                index_d = i;
            end
        end
    end
    % Next, find smallest time to next pulse firing
    [firingtime, index_f] = min((threshold - theta)./omega);
    % Find overall smallest time between these two values
    [timestep, index] = min([receivetime, firingtime]);
    
    %%Advance simulation to next event
    % Updated phases to next event
    theta = theta + timestep*omega;
    % Decrease delays by timestep
    for i = 1:1:length(delay)
        if isempty(delay{i}) == false % If oscillator i has a delay timer
            delay{i} = delay{i} - timestep;
        end
    end
    % Store data for plotting evolution of phases for the episode
    k = k + 1; % Increment data point number
    time(k) = time(k-1) + timestep; % Current time
    state(k,:) = theta; % Current state of network
    c(k) = ContainingArc(state(k,:), threshold); % Current Containing Arc of network
    for i = 1:N
        % Current estimate of containing arc for each oscillator
        c_est(k,i) = ContainingArc([delta_matrix{i}, 0],threshold);
    end
    
    %%Update simulation state based on the event
    d_theta = zeros(size(theta)); % Change in phase based on the event
    if index == 1 % Next event is oscillator index_d responding to received pulse
        % If theta(index_d) < D % If phase is inside refractory period
            % Then skip this; no response for this event; remove delay timer
        % end
        % Else, with probability P_drop (i.e., dropped pulse), skip this;
        % no response for this event; remove delay timer
        
        % Determine action based on episode's PRC
        Action = GetAction(theta(index_d),PRC(:,index_d),State_positions);
        d_theta(index_d) = l * Action; % Desired change in phase for oscillator index_d
        % Remove timer
        delay{index_d}(j) = [];
        
        % Oscillator index_d updates phase difference estimates of other oscillators
        % m is the number of received pulses since its own firing
        m(index_d) = m(index_d) + 1; % Increment received pulse counter
        delta_matrix{index_d}(m(index_d)) = mod(-(theta(index_d) + d_theta(index_d)),threshold); % Update for oscillator that fired
        for i = 1:length(delta_matrix{index_d})
            %disp("Testing")
            if i == m(index_d) % delta_matrix{index_d}(i) corresponds to oscillator that fired a pulse
                %delta_matrix{index_d}(i) = mod(-(theta(index_d) + d_theta(index_d)),threshold); % Update for oscillator that fired
            else %if i ~= m(index_d) % For all other oscillators; assume pulse was received at the same time
                % Change to phase difference estimate
                % update based on current PRC
                %stuff = l * GetAction(mod(theta(index_d)+delta_matrix{index_d}(i), threshold), PRC, State_positions) - d_theta(index_d);
                % update based on stale estimate
                stuff = 0 - d_theta(index_d);
                delta_matrix{index_d}(i) = mod(delta_matrix{index_d}(i) + stuff, 2*pi);
            end
        end
        
        % Record S-A pairs (with probability) for learning
        SA_pairs(index_d, :, k_SAR(index_d)) = [theta(index_d), Action];
        
    elseif index == 2 % Next event is pulse firing by oscillator index_f
        % Determine action based on episode's PRC
        %Action = GetAction(theta(index_f)-threshold,PRC,State_positions);
        %d_theta(index_f) = -threshold + l * Action; % Recent phase of firing oscillator
        d_theta(index_f) = -threshold;
        % Add delays to all connected oscillators
        for i = 1:1:N
            if A(index_f,i) == 1 % If oscillator i is connected to oscillator index_f...
                % With probability 1... (no pulse drops)
                % Add delay timer to oscillator i for next received pulse
                %delay{i} = [delay{i}, pulse_delay(index_f,i)]; % Fixed delay
                delay{i} = [delay{i}, d_c + (d_v)*rand(1)]; % Random delay
            end
        end
        
        if m(index_f) < length(delta_matrix{index_f}) % If pulses received is less than state differences that oscilator index_f is tracking...
            % Remove excess state differences
            delta_matrix{index_f}(m(index_f)+1:end) = [];
        end % Otherwise, leave delta_matrix alone
        m(index_f) = 0; % reset received pulse counter
        
        % Oscillator index_f updates phase difference estimates of all oscillators
        %for i = 1:length(delta_matrix{index_f})
            % Change to phase difference estimate
            %stuff = l * GetAction(mod(delta_matrix{index_f}(i), threshold), PRC, State_positions);
            %delta_matrix{index_f}(i) = delta_matrix{index_f}(i) + stuff;
        %end
        
        % Record S-A pairs (with probability) for learning
        %SA_pairs(index_f, :, k_SAR(index_f)) = [theta(index_f)-threshold, Action];
    end
    
    % Updated phases based on previous event
    theta = mod(theta + d_theta, threshold); % Mod should be unnecessary
    % Store data
    k = k + 1;
    time(k) = time(k-1);
    state(k,:) = theta;
    c(k) = ContainingArc(state(k,:), threshold);
    for i = 1:N
        c_est(k,i) = ContainingArc([0, delta_matrix{i}],threshold);
    end
    
    energy_option = 1;
    w_Lambda = N/(N-1); % Weighting of containing arc in reward value
    w_Energy = (1/l); % Weighting of action energy in reward value
    % Determine reward; determined at firing instance of oscillator
    if index == 1 % If last event was a pulse receiving...
        % Calculate energy used to take Action
        %w_Lambda = ((in_degree(index_d)+1)/(in_degree(index_d))); % Weighting of containing arc in reward value
        r_Energy = GetRewardEnergy(energy_option,d_theta(index_d),threshold);
        
        % Compute reward; = (Change in containing arc from previous firing) 
        % - (enery expended since previous firing)
        R(k_SAR(index_d),index_d) = w_Lambda*(r_Lambda(index_d) - c_est(k,index_d)) - w_Energy*r_Energy;
        %R(k_SAR(index_d),index_d) = (N/(N-1))*(r_Lambda(index_d) - c(k)) - r_Energy;
        r_Lambda(index_d) = c_est(k,index_d); % Save current current containing arc for next firing instance
        %r_Lambda(index_d) = c(k); % Save current current containing arc for next firing instance
        k_SAR(index_d) = k_SAR(index_d) + 1; % Increment cycle counter
    %elseif index == 2
        % Calculate energy used to take Action
        %r_Energy = GetRewardEnergy(energy_option,d_theta(index_f)+threshold,threshold);
        
        % Compute reward; = (Change in containing arc from previous firing) 
        % - (enery expended since previous firing)
        %R(k_SAR(index_f),index_f) = (N/(N-1))*(r_Lambda(index_f) - c_est(k,index_f)) - r_Energy;
        %R(k_SAR(index_f),index_f) = (N/(N-1))*(r_Lambda(index_f) - c(k)) - r_Energy;
        %r_Lambda(index_f) = c_est(k,index_f); % Save current current containing arc for next firing instance
        %r_Lambda(index_f) = c(k); % Save current current containing arc for next firing instance
        %k_SAR(index_f) = k_SAR(index_f) + 1; % Increment cycle counter
    end
end % while time(k) < T_eps
% Perhaps add condition if PCO simulation stalls;
end

% Calculate reward penalty due to oscillator action
function [energy] = GetRewardEnergy(option,Action,threshold)
% Calculate energy used to take Action
switch option
    case 1 % Quadratic
        energy = (Action^2)/threshold; %(Action^2)/threshold;
    case 2 % Logarithmic
        energy = (threshold)*log1p((abs(Action)/threshold));
    case 3 % Average of Quadratic and Logarthmic
        energy = 0.5*((Action^2)/threshold + threshold*log1p((abs(Action)/threshold)));
    otherwise
        % Others
        %energy = sqrt(abs(Action));
        energy = abs(Action);
        %energy = (1/2)*(sqrt(abs(Action)) + Action^2);
        %energy = (Action^2) + 1;
        %energy = (abs(Action)) + 1;
end
end

% Update the State-Action values based on the episode results
% Parameters:
%   Q_old = matrix; previous Q-values for state-action pairs
%   SA_pairs = matrix; State-actions used in previous episode
%   R = matrix; Reward for each oscillator based on actions in previous episode
%   alpha_SAR = float; SARSA learning rate (0,1)
%   gamma_SAR = float; SARSA discount rate (0,1)
% Returns:
%   Q_new = matrix; updated Q-values for state-action pairs
function [Q_new] = QUpdateSARSA(Q_old, SA_pairs, R, alpha_SAR, gamma_SAR, State_positions, Action_space, PRC, threshold)
[~, num_Action] = size(Action_space);
deltaQ = zeros(size(Q_old)); % Change to Q value for each state-action pair
for n = 1:length(SA_pairs(:,1,1)) % For each oscillator in the network
    for k = 1:length(SA_pairs(n,1,:))-1 % For each state-action pair
        if SA_pairs(n,:,k) == zeros(size(SA_pairs(n,:,k)))
            continue % No pulses were received before firing
        end
        % Else, Determine SARSA update to Q-values for S-A pair
        R1 = R(k,n); % determine reward for these actions
        % Get low and high state parameter indices for current S-A pair (SA_pairs(i,1,k))
        [state1_low, S1_low] = max(State_positions(State_positions < SA_pairs(n,1,k)));
        % Catch condition if SA_pairs(i,1,k) == 0
        if isempty(state1_low) % Indicates that SA_pairs(i,1,k) == 0
            S1_low = 1;
            state1_low = State_positions(S1_low);
        end
        S1_high = S1_low + 1;
        state1_high = State_positions(S1_high);
        % Get low and high action indices for current S-A pair (SA_pairs(i,2,k))
        % General procedure:
        %A1_low_possible = find(Action_space(S1_low,:) == PRC(S1_low,n));
        %A1_high_possible = find(Action_space(S1_high,:) == PRC(S1_high,n));
        %A1_low = A1_low_possible(randsample(length(A1_low_possible),1));
        %A1_high = A1_high_possible(randsample(length(A1_high_possible),1));
        % Direct calculation based on discretation of action space
        A1_low = round(((PRC(S1_low,n) + state1_low)/threshold) * (num_Action-1)) + 1;
        A1_high = round(((PRC(S1_high,n) + state1_high)/threshold) * (num_Action-1)) + 1;
        % "Closeness" to low and high state parameters
        low_prob1 = (state1_high - SA_pairs(n,1,k))/(state1_high - state1_low);
        high_prob1 = 1 - low_prob1;
        
        % Get low and high state parameter indices for current S-A pair (SA_pairs(i,1,k+1))
        [state2_low, S2_low] = max(State_positions(State_positions < SA_pairs(n,1,k+1)));
        % Catch condition if SA_pairs(i,1,k+1) == 0
        if isempty(state2_low) % Indicates that SA_pairs(i,1,k+1) == 0
            S2_low = 1;
            state2_low = State_positions(S2_low);
        end
        S2_high = S2_low + 1;
        state2_high = State_positions(S2_high);
        % Get low and high action indices for current S-A pair (SA_pairs(i,2,k+1))
        % General procedure
        %A2_low_possible = find(Action_space(S2_low,:) == PRC(S2_low,n));
        %A2_high_possible = find(Action_space(S2_high,:) == PRC(S2_high,n));
        %A2_low = A2_low_possible(randsample(length(A2_low_possible),1));
        %A2_high = A2_high_possible(randsample(length(A2_high_possible),1));
        % Direct calculation based on discretation of action space
        A2_low = round(((PRC(S2_low,n) + state2_low)/threshold) * (num_Action-1)) + 1;
        A2_high = round(((PRC(S2_high,n) + state2_high)/threshold) * (num_Action-1)) + 1;
        % "Closeness" to low and high state parameters
        low_prob2 = (state2_high - SA_pairs(n,1,k+1))/(state2_high - state2_low);
        high_prob2 = 1 - low_prob2;
        
        % Determine value of [bracket] term (based on combined high/low state-action values
        bracket = low_prob2*Q_old(S2_low,A2_low,n) + high_prob2*Q_old(S2_high,A2_high,n);
        %{
        for m = 1:length(SA_pairs(i,1,:,k))
            if SA_pairs(i,2,:,k+1) == 0
                continue % Does not contribute to bracket term
            end
            %[S2,A2,Sprob2] = SA_pairs(i,:,m,k+1);
            S2 = SA_pairs(i,1,m,k+1);
            A2 = SA_pairs(i,2,m,k+1); % A2 = PRC(S2)
            Sprob2 = SA_pairs(i,3,m,k+1);
            bracket = bracket + Sprob2*Q_old(S2,A2);
        end
        %}
        
        update_low = low_prob1*alpha_SAR*(R1 + gamma_SAR*bracket - Q_old(S1_low,A1_low,n));
        deltaQ(S1_low,A1_low,n) = deltaQ(S1_low,A1_low,n) + update_low;
        update_high = high_prob1*alpha_SAR*(R1 + gamma_SAR*bracket - Q_old(S1_high,A1_high,n));
        deltaQ(S1_high,A1_high,n) = deltaQ(S1_high,A1_high,n) + update_high;
        %{
        for m = 1:length(SA_pairs(i,1,:,k))
            % What does this condition do exactly?
            if SA_pairs(i,:,m,k) == zeros(size(SA_pairs(i,:,m,k)))
                continue % Not a pair - Extra padding in SA_pair array
            end
            %[S1,A1,Sprob1] = SA_pairs(i,:,m,k); % Grab SA values
            S1 = SA_pairs(i,1,m,k);
            A1 = SA_pairs(i,2,m,k); % A1 = PRC(S1)
            Sprob1 = SA_pairs(i,3,m,k);
            
            % SARSA Update to Q-values for S-A pair
            update = alpha_SAR*Sprob1*(R1 + gamma_SAR*bracket - Q_old(S1,A1));
            deltaQ(S1,A1) = deltaQ(S1,A1) + update;
        end
        %}
    end % End for cycle k, move to next cycle
end % End of oscillator i, move to next oscillator

Q_new = Q_old + deltaQ; % Perform update to Q-values from episode
% Q(S,A) = Q(S,A) + prob*alpha*(R + gamma*(bracket) - Q(S,A))
% SA_pairs(i,:,k) goes with R(k) and SA_pairs(i,:,k+1)
% If SA_pairs(i,:,k) == [0,0,0], skip pair.
end

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
    case 'star' % Star
        A = MakeTopology(N, 2:N);
        A(2:N,2:N) = zeros(N-1);
    case 'ERGrand' % Erdős–Rényi–Gilbert random network
        A = rand(N,N) + eye(N);
        if symmetric % birectional graph
            A = A - tril(A,-1); % remove lower triangular components
            A = A+A';
        end
        A = double(A<vector); % vector is probability of making an edge (0,1]
    case 'alldel' % special matrix; nodes have every possible indegree
        if symmetric
            halfi = triu(rot90(eye(N)));
            A = rot90(triu(ones(N))-halfi);
        else
            A = triu(ones(N),1);
            %A(2,1) = 1;
        end
    case 'custom' % Custom
        A = MakeTopology(N, vector, symmetric); % Custom
        % Include other modifications to A here.
    otherwise % Nothing
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

function DisplayPercentageComplete(iter,total,time,step)
percent = 100*iter/total;
for i = 0:step:100
    if percent == i
        disp([num2str(i),'% Complete in ',num2str(time),' seconds'])
        break % We done
    end
end
end

function [BigQ] = ScaleQ(Q)
[Pq, Aq] = size(Q);
P = Pq - 1;
P_prime = 2*P + 1;
A = Aq - 1;
A_prime = 2*A + 1;
BigQ = zeros(P_prime,A_prime);
for i = 1:1:Pq
    for j = 1:1:Aq
        BigQ((2*i)-1,(2*j)-1) = Q(i,j);
    end
end
for i = 1:1:P
    for j = 1:1:Aq
        BigQ((2*i),(2*j)-1) = 0.5 * (Q(i,j) + Q(i+1,j));
    end
end
for i = 1:1:Pq
    for j = 1:1:A
        BigQ((2*i)-1,(2*j)) = 0.5 * (Q(i,j) + Q(i,j+1));
    end
end
for i = 1:1:P
    for j = 1:1:A
        BigQ((2*i),(2*j)) = 0.25 * (Q(i,j) + Q(i+1,j) + Q(i,j+1) + Q(i+1,j+1));
    end
end
end