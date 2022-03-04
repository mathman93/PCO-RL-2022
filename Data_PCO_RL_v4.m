% Analytical Data from simulations based on PCO_RL_v4.m
% Optimal PRFs with heterogeneous state-action values
% Use in conjunction with PCO_RL_v4.m and GeneralPRFCheck.m

%% Save Table
writetable(RLPRC_data,'RLPRF_data_rand.txt', 'Delimiter', ',')

%% Create Table of Data
disp("Do you want to do this?") % Catch for accidently deleting data
pause
table_vartypes = {'int32', 'int32', 'string', 'double', 'int32', 'double', 'double', 'double', 'double', 'double'};
table_varnames = {'N', 'OscNum', 'Topology', 'Coupling', 'InDegree', 'FreqVar', 'DelayConst', 'DelayVar', 'c1', 'c2'};
RLPRC_data = table('Size', [0,10], 'VariableTypes', table_vartypes, 'VariableNames', table_varnames);
% This is an empty table

%% Load Table
RLPRC_data = readtable('RLPRF_data_rand.txt', 'Delimiter', ',', 'Format', '%u%u%s%f%u%f%f%f%f%f');
disp(RLPRC_data)

%% Display Entire Table
disp(RLPRC_data)

%% Sort table
RLPRC_data = sortrows(RLPRC_data, [5,1,3,2]);
% Sort by indegree, N, topology, then oscillator #.
disp(RLPRC_data)

%% Add Data to Table
h = height(RLPRC_data);
for i = 1:N
    RLPRC_data(h+i,:) = {N, i, char(A_option), l, in_degree(i), freq_var, delay_const, delay_var, c1(i)*(threshold/P), c2(i)*(threshold/P)};
end
disp(RLPRC_data((h+1):(h+N),:))
disp(height(RLPRC_data))

%% Remove Data from Table
rows = 0; % Range of rows to remove
RLPRC_data(rows,:) = [];
disp(RLPRC_data)

%% Data Analysis

%Non-linear best fit (inverse w/ offset) line for average values
modelfunc = @(b,x)((2*pi)*(((1+b(1))*(0.5-b(2))./(x+b(1)))+b(2)));
%Non-linear best fit (exponential w/ offset) line for average values
modelfunc1 = @(b,x)((2*pi)*((0.5-b(2))*exp(-b(1)*(x-1))+b(2)));
%modelfunc2 = @(b,x)(2*pi-((pi-b(1))*exp(-b(2)*(x-1))+b(1)));
%opts = statset('nlinfit');
%opts.RobustWgtFun = 'bisquare';
%beta0 = [0,pi];

% Combine data based on topology (symmetric networks)
%{
labels = sortrows(unique(RLPRC_data(:,[1,3,5])),[3,1,2]); % Get unique topologies found in data set
% Allocate space
h_labels = height(labels);
lej = strings(h_labels,1);
delt = zeros(h_labels,1);
c1_avg = zeros(h_labels,1);
c1_max = zeros(h_labels,1);
c1_min = zeros(h_labels,1);
c2_avg = zeros(h_labels,1);
c2_max = zeros(h_labels,1);
c2_min = zeros(h_labels,1);
for i = 1:h_labels % For each unique topology...
    % Separate out label info
    N_label = labels.N(i);
    topology_label = labels.Topology(i);
    delt(i) = labels.InDegree(i);
    % Find data for this specific label
    mini_data = RLPRC_data.N == N_label & RLPRC_data.Topology == string(topology_label);
    lej(i) = LegendString(N_label, topology_label); % Legend Info
    %lej(i) = "N = " + string(N_label) + "; " + string(topology_label); % Legend Info
    %disp(RLPRC_data(mini_data,:))
    % PRF parameter values for model fitting and stats
    c1_vals = RLPRC_data.c1(mini_data);
    c2_vals = RLPRC_data.c2(mini_data);
    c1_avg(i) = mean(c1_vals);
    c1_max(i) = max(c1_vals);
    c1_min(i) = min(c1_vals);
    c2_avg(i) = mean(c2_vals);
    c2_max(i) = max(c2_vals);
    c2_min(i) = min(c2_vals);
end
%}
% Use data point individually (random networks)
%
the_data = RLPRC_data.Topology == "ERGrand";
labels = sortrows(unique(RLPRC_data(:,[5])),[1]); % Get unique topologies found in data set
% Allocate space
h_labels = height(labels);
%lej = strings(h_labels,1);
delt = zeros(h_labels,1);
c1_avg = zeros(h_labels,1);
c1_max = zeros(h_labels,1);
c1_min = zeros(h_labels,1);
c2_avg = zeros(h_labels,1);
c2_max = zeros(h_labels,1);
c2_min = zeros(h_labels,1);
for i = 1:h_labels % For each unique indegree value...
    % Separate out label info
    %N_label = labels.N(i);
    %topology_label = labels.Topology(i);
    delt(i) = labels.InDegree(i);
    % Find data for this specific label
    mini_data = RLPRC_data.InDegree == delt(i); %& RLPRC_data.Topology == "ERGrand";
    %lej(i) = LegendString(N_label, topology_label); % Legend Info
    %lej(i) = "N = " + string(N_label) + "; " + string(topology_label); % Legend Info
    %disp(RLPRC_data(mini_data,:))
    % PRF parameter values for model fitting and stats
    c1_vals = RLPRC_data.c1(mini_data);
    c2_vals = RLPRC_data.c2(mini_data);
    c1_avg(i) = mean(c1_vals);
    c1_max(i) = max(c1_vals);
    c1_min(i) = min(c1_vals);
    c2_avg(i) = mean(c2_vals);
    c2_max(i) = max(c2_vals);
    c2_min(i) = min(c2_vals);
end
%}
% Determine best fit model parameters
beta1 = nlinfit(double(delt),(c1_avg),modelfunc, [3,0.1]);
beta2 = nlinfit(double(delt),(c2_avg),modelfunc, [3,0.9]);
beta3 = nlinfit(double(delt),(c1_avg),modelfunc1, [0.5,0.1]);
beta4 = nlinfit(double(delt),(c2_avg),modelfunc1, [0.5,0.9]);
disp([beta1, beta2])
disp([beta3, beta4])

% Define points for best fit models
x_fit = 1:0.1:100;
y1_fit = modelfunc(beta1,x_fit);
y2_fit = modelfunc(beta2,x_fit);
y3_fit = modelfunc1(beta3,x_fit);
y4_fit = modelfunc1(beta4,x_fit);

%{
fig2h = figure(2);
clf
set(fig2h, 'Position', [800, 60, 450, 480])
markers = 'xo<*s>';
sub1h = subplot(2,1,1);
%set(sub1h, 'TickLabelInterpreter', 'latex')
set(sub1h, 'Position', [0.12,0.59, 0.81, 0.40])
hold on
c2_neg = c2_avg - c2_min;
c2_pos = c2_max - c2_avg;
errorbar(delt, c2_avg, c2_neg, c2_pos, 'bx', 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)

%plot(x_fit,y2_fit, 'r--', 'Linewidth', 1.5)
plot(x_fit,y4_fit, 'k:', 'Linewidth', 1.5)
fitline2 = annotation('textarrow',[0.42,0.45],[0.94,0.86],...
    'String','$c_{2} \approx (-0.51\pi)e^{-0.14(\delta^{-}-1)} + 1.51\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on
xlabelh = xlabel('');
ylabelh = ylabel('$c_{2}$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14, 'Rotation', 0)
axis([0.1,23,0.45*threshold,0.80*threshold])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

sub2h = subplot(2,1,2);
%set(sub2h, 'TickLabelInterpreter', 'latex')
set(sub2h, 'Position', [0.12,0.08, 0.81, 0.46])
hold on   
c1_neg = c1_avg - c1_min;
c1_pos = c1_max - c1_avg;
errorbar(delt, c1_avg, c1_neg, c1_pos, 'bx', 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)

%plot(x_fit,y1_fit, 'r--', 'Linewidth', 1.5)
plot(x_fit,y3_fit, 'k:', 'Linewidth', 1.5)
fitline1 = annotation('textarrow',[0.41,0.44],[0.14,0.19],...
    'String','$c_{1} \approx (0.85\pi)e^{-0.23(\delta^{-}-1)} + 0.15\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on
xlabelh = xlabel('Network Indegree $(\delta^{-})$', 'interpreter', 'latex');
ylabelh = ylabel('$c_{1}$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14, 'Rotation', 0)
axis([0.1,23,0,0.54*threshold])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

%}

% Figure 2 alternate (combined)
%
c1_neg = c1_avg - c1_min;
c1_pos = c1_max - c1_avg;
c2_neg = c2_avg - c2_min;
c2_pos = c2_max - c2_avg;

fig2h = figure(2);
clf
set(fig2h, 'Position', [800, 60, 420, 320])
hold on
errorbar(delt, c2_avg, c2_neg, c2_pos, 'x', 'Color', [0.6,0.15,0.9], 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)
errorbar(delt, c1_avg, c1_neg, c1_pos, 'x', 'Color', [0.95,0.4,0], 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)

%plot(x_fit,y2_fit, 'r--', 'Linewidth', 1.5)
plot(x_fit,y4_fit, 'k-.', 'Linewidth', 1.5)
fitline2 = annotation('textarrow',[0.43,0.45],[0.865,0.80],...
    'String','$c_{2} \approx (-0.51\pi)e^{-0.14(\delta^{-}-1)} + 1.51\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
%plot(x_fit,y1_fit, 'r:', 'Linewidth', 1.5)
plot(x_fit,y3_fit, 'k:', 'Linewidth', 1.5)
fitline1 = annotation('textarrow',[0.42,0.44],[0.20,0.265],...
    'String','$c_{1} \approx (0.81\pi)e^{-0.29(\delta^{-}-1)} + 0.19\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on

xlabelh = xlabel('Oscillator Indegree $(\delta^{-})$', 'interpreter', 'latex');
ylabelh = ylabel('Parameter Values', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14)%, 'Rotation', 0)
axis([0.2,21,0,0.8*threshold])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

legendh = legend('$c_2$', '$c_1$');%, 'interpeter', 'latex')
set(legendh,'Fontname','Times New Roman', 'Fontsize',11, 'Interpreter', 'latex', 'Position', [0.72, 0.42, 0.16, 0.14])

%}

%% OLD Data Analysis (for RLPRF_data.txt)

%Non-linear best fit (inverse w/ offset) line for average values
modelfunc = @(b,x)((2*pi)*(((1+b(1))*(0.5-b(2))./(x+b(1)))+b(2)));
%Non-linear best fit (exponential w/ offset) line for average values
modelfunc1 = @(b,x)((2*pi)*((0.5-b(2))*exp(-b(1)*(x-1))+b(2)));
%modelfunc2 = @(b,x)(2*pi-((pi-b(1))*exp(-b(2)*(x-1))+b(1)));
%opts = statset('nlinfit');
%opts.RobustWgtFun = 'bisquare';
%beta0 = [0,pi];

% Combine data based on topology (symmetric networks)
labels = sortrows(unique(RLPRC_data(:,[1,3,5])),[3,1,2]); % Get unique topologies found in data set
% Allocate space
h_labels = height(labels);
lej = strings(h_labels,1);
delt = zeros(h_labels,1);
c1_avg = zeros(h_labels,1);
c1_max = zeros(h_labels,1);
c1_min = zeros(h_labels,1);
c2_avg = zeros(h_labels,1);
c2_max = zeros(h_labels,1);
c2_min = zeros(h_labels,1);
for i = 1:h_labels % For each unique topology...
    % Separate out label info
    N_label = labels.N(i);
    topology_label = labels.Topology(i);
    delt(i) = labels.InDegree(i);
    % Find data for this specific label
    mini_data = RLPRC_data.N == N_label & RLPRC_data.Topology == string(topology_label);
    lej(i) = LegendString(N_label, topology_label); % Legend Info
    %lej(i) = "N = " + string(N_label) + "; " + string(topology_label); % Legend Info
    %disp(RLPRC_data(mini_data,:))
    % PRF parameter values for model fitting and stats
    c1_vals = RLPRC_data.c1(mini_data);
    c2_vals = RLPRC_data.c2(mini_data);
    c1_avg(i) = mean(c1_vals);
    c1_max(i) = max(c1_vals);
    c1_min(i) = min(c1_vals);
    c2_avg(i) = mean(c2_vals);
    c2_max(i) = max(c2_vals);
    c2_min(i) = min(c2_vals);
end

% Determine best fit model parameters
beta1 = nlinfit(double(delt),(c1_avg),modelfunc, [3,0.1]);
beta2 = nlinfit(double(delt),(c2_avg),modelfunc, [3,0.9]);
beta3 = nlinfit(double(delt),(c1_avg),modelfunc1, [0.5,0.1]);
beta4 = nlinfit(double(delt),(c2_avg),modelfunc1, [0.5,0.9]);
disp([beta1, beta2])
disp([beta3, beta4])

% Define points for best fit models
x_fit = 1:0.1:100;
y1_fit = modelfunc(beta1,x_fit);
y2_fit = modelfunc(beta2,x_fit);
y3_fit = modelfunc1(beta3,x_fit);
y4_fit = modelfunc1(beta4,x_fit);

%{
dataset1 = RLPRC_data.Coupling == 0.99 & RLPRC_data.FreqVar == 0 & RLPRC_data.DelayConst == 0 & RLPRC_data.DelayVar == 0;
disp(sum(dataset1))
delta = RLPRC_data.InDegree(dataset1);% & RLPRC_data.N < 10);
c1_val = RLPRC_data.c1(dataset1);% & RLPRC_data.N < 10);
c2_val = RLPRC_data.c2(dataset1);% & RLPRC_data.N < 10);
beta1 = nlinfit(double(delta),(c1_val),modelfunc, [3,0.1]);
beta2 = nlinfit(double(delta),(c2_val),modelfunc, [3,0.9]);
beta3 = nlinfit(double(delta),(c1_val),modelfunc1, [0.5,0.1]);
beta4 = nlinfit(double(delta),(c2_val),modelfunc1, [0.5,0.9]);
disp([beta1, beta2])
disp([beta3, beta4])

x_fit = 1:0.1:100;
y1_fit = modelfunc(beta1,x_fit);
y2_fit = modelfunc(beta2,x_fit);
y3_fit = modelfunc1(beta3,x_fit);
y4_fit = modelfunc1(beta4,x_fit);

%RLPRC_data.N == 3
x_val = RLPRC_data.InDegree(dataset1);% & RLPRC_data.N < 10);
x_unique = unique(x_val);
c1_avg = zeros(size(x_unique));
c1_max = zeros(size(x_unique));
c1_min = zeros(size(x_unique));
c2_avg = zeros(size(x_unique));
c2_max = zeros(size(x_unique));
c2_min = zeros(size(x_unique));
for i = 1:length(x_unique)
    n = x_unique(i);
    c1_data = RLPRC_data.c1(RLPRC_data.Coupling == 0.99 & RLPRC_data.InDegree == n);% & RLPRC_data.N < 10);
    c1_avg(i) = mean(c1_data);
    c1_max(i) = max(c1_data);
    c1_min(i) = min(c1_data);
    c2_data = RLPRC_data.c2(RLPRC_data.Coupling == 0.99 & RLPRC_data.InDegree == n);% & RLPRC_data.N < 10);
    c2_avg(i) = mean(c2_data);
    c2_max(i) = max(c2_data);
    c2_min(i) = min(c2_data);
end
%}

%
fig2h = figure(2);
clf
set(fig2h, 'Position', [800, 60, 450, 480])
markers = 'xo<*s>';
sub1h = subplot(2,1,1);
%set(sub1h, 'TickLabelInterpreter', 'latex')
set(sub1h, 'Position', [0.12,0.59, 0.81, 0.40])
hold on
%plot(RLPRC_data.InDegree(RLPRC_data.CouplingStrength == 0.99),...
%     RLPRC_data.c2(RLPRC_data.CouplingStrength == 0.99), 'kx')
for i = 1:h_labels
    c2_neg = c2_avg(i) - c2_min(i);
    c2_pos = c2_max(i) - c2_avg(i);
    errorbar(delt(i), c2_avg(i), c2_neg, c2_pos,...
        markers(mod(i-1,length(markers))+1), 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)
end
%plot(x_fit,y2_fit, 'r--', 'Linewidth', 1.5)
plot(x_fit,y4_fit, 'k:', 'Linewidth', 1.5)
fitline2 = annotation('textarrow',[0.42,0.45],[0.94,0.86],...
    'String','$c_{2} \approx (-0.53\pi)e^{-0.13(\delta^{-}-1)} + 1.53\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on
xlabelh = xlabel('');
ylabelh = ylabel('$c_{2}$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14, 'Rotation', 0)
axis([0.1,23,0.45*threshold,0.80*threshold])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

sub2h = subplot(2,1,2);
%set(sub2h, 'TickLabelInterpreter', 'latex')
set(sub2h, 'Position', [0.12,0.08, 0.81, 0.46])
hold on
%plot(RLPRC_data.InDegree(RLPRC_data.CouplingStrength == 0.99),...
%     RLPRC_data.c1(RLPRC_data.CouplingStrength == 0.99), 'kx')   
for i = 1:h_labels
    c1_neg = c1_avg(i) - c1_min(i);
    c1_pos = c1_max(i) - c1_avg(i);
    errorbar(delt(i), c1_avg(i), c1_neg, c1_pos,...
        markers(mod(i-1,length(markers))+1), 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)
end
%plot(x_fit,y1_fit, 'r--', 'Linewidth', 1.5)
plot(x_fit,y3_fit, 'k:', 'Linewidth', 1.5)
fitline1 = annotation('textarrow',[0.41,0.44],[0.14,0.19],...
    'String','$c_{1} \approx (0.85\pi)e^{-0.23(\delta^{-}-1)} + 0.15\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on
xlabelh = xlabel('Network Indegree $(\delta^{-})$', 'interpreter', 'latex');
ylabelh = ylabel('$c_{1}$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14, 'Rotation', 0)
axis([0.1,23,0,0.54*threshold])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

legendh = legend(lej);
set(legendh,'Fontname','Times New Roman', 'Fontsize',8, 'Interpreter', 'latex', 'Position', [0.72, 0.33, 0.23, 0.40])

%legendh = legend('$\delta^{-}=1$', '$\delta^{-}=2$', '$\delta^{-}=3$', '$\delta^{-}=4$',...
%    '$\delta^{-}=5$', '$\delta^{-}=6$', '$\delta^{-} = 10$', '$\delta^{-} = 15$');%, '$\delta^{-} = 20$');
%set(legendh,'Fontname','Times New Roman', 'Fontsize',10, 'Interpreter', 'latex', 'Position', [0.72, 0.34, 0.23, 0.40])
%}

%% Prior Analysis
%{
% Determine best fit (exponential) line for average values
A_m = delta-1;
y1 = log(pi) - log(c1_avg*2*pi);
y2 = log(pi) - log(2*pi - c2_avg*2*pi);
s1 = (A_m'*y1)/(A_m'*A_m); %Normal Equation
s2 = (A_m'*y2)/(A_m'*A_m);
x_val = 1:.1:20;
% c1 = pi*exp(-s1*(delta-1)); c2 = 2*pi - pi*exp(-s2*(delta-1));

% Determine best fit (inverse proportional) line for average values
A_m2 = log(delta);
y1 = log(pi) - log(c1_avg*2*pi);
y2 = log(pi) - log(2*pi - c2_avg*2*pi);
r1 = (A_m2'*y1)/(A_m2'*A_m2); %Normal Equation
r2 = (A_m2'*y2)/(A_m2'*A_m2);
%x_val = 0:0.1:19;
% c1 = pi*(delta+1)^(-r1); c2 = 2*pi - pi*(delta+1)^(-r2));

%Non-linear best fit (exponential w/ offset) line for average values
modelfunc1 = @(b,x)((pi-b(1))*exp(-b(2)*(x-1))+b(1));
modelfunc2 = @(b,x)(2*pi-((pi-b(1))*exp(-b(2)*(x-1))+b(1)));
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
beta0 = [0,0.5];
beta1 = nlinfit((delta),(c1_avg*2*pi),modelfunc1, beta0);
beta2 = nlinfit((delta),(c2_avg*2*pi),modelfunc2, beta0);

%{
fig2h = figure(2);
clf
set(fig2h, 'Position', [850, 60, 420, 500])

markers = 'xo<*s>';

sub1h = subplot(2,1,1);
%set(sub1h, 'TickLabelInterpreter', 'latex')
set(sub1h, 'Position', [0.13,0.56, 0.81, 0.43])
hold on
for i = 1:length(delta)
    k = delta_index(i);
    errorbar(delta(k), 2*pi*c2_avg(k), 2*pi*(c2_avg(k)-c2_min(k)), 2*pi*(c2_max(k)-c2_avg(k)),...
        markers(mod(i-1,length(markers))+1), 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)
end
%plot(x_val+1, 2*pi - pi*exp(-s2*x_val), 'r--', 'Linewidth', 1.5)
%plot(x_val+1, 2*pi - pi*(x_val+1).^(-r2), 'b-.', 'Linewidth', 1.5)
plot(x_val, modelfunc2(beta2,x_val), 'k:', 'Linewidth', 1.5)
fitline2 = annotation('textarrow',[0.60,0.64],[0.92,0.84],...
    'String','$c_{2} \approx 1.60\pi - (0.60\pi)e^{-0.26(\delta^{-}-1)}$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on
xlabelh = xlabel('');
ylabelh = ylabel('$c_{2}$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14, 'Rotation', 0)
axis([0.5,10.5,3.04,5.54])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

sub2h = subplot(2,1,2);
%set(sub2h, 'TickLabelInterpreter', 'latex')
set(sub2h, 'Position', [0.13,0.08, 0.81, 0.43])
hold on
for i = 1:length(delta)
    k = delta_index(i);
    errorbar(delta(k),2*pi*c1_avg(k), 2*pi*(c1_avg(k)-c1_min(k)), 2*pi*(c1_max(k)-c1_avg(k)),...
        markers(mod(i-1,length(markers))+1), 'Linewidth', 1.0, 'CapSize', 8, 'MarkerSize', 6)
end
%plot(x_val+1, pi*exp(-s1*x_val), 'r--', 'Linewidth', 1.5)
%plot(x_val+1, pi*(x_val+1).^(-r1), 'b-.', 'Linewidth', 1.5)
plot(x_val, modelfunc1(beta1,x_val), 'k:', 'Linewidth', 1.5)
fitline1 = annotation('textarrow',[0.66,0.71],[0.125,0.145],...
    'String','$c_{1} \approx (0.72\pi)e^{-0.33(\delta^{-}-1)} + 0.28\pi$',...
    'FontName','Times New Roman', 'Fontsize',12, 'headStyle','vback3', 'lineStyle','-', 'HeadLength',7, 'HeadWidth',7, 'Interpreter','latex');
hold off
grid on
xlabelh = xlabel('Network Indegree $(\delta^{-})$', 'interpreter', 'latex');
ylabelh = ylabel('$c_{1}$', 'interpreter', 'latex');
set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',14, 'Rotation', 0)
axis([0.5,10.5,0.74,3.24])
yticks([0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

legendh = legend('All-to-All (N=2)',...
    'All-to-All (N=3)', 'Ring (N=4)', 'Ring (N=5)', 'Ring (N=6)', 'Ring (N=7)',...
    'All-to-All (N=4)', 'Near2-Far1 (N=6)', 'Far3 (N=6)',...
    'All-to-All (N=5)', 'Near4 (N=6)', 'Near4 (N=7)',...
    'All-to-All (N=6)', 'All-to-All (N=7)', 'All-to-All (N=11)',...
    'Best Fit Trendline');
set(legendh,'Fontname','Times New Roman', 'Fontsize',10, 'Position', [0.62, 0.24, 0.36, 0.54])
%}

%}
%% Graph of Learned Function (for various network indegrees)
%
PRF_gx = zeros(7,5);
PRF_gy = zeros(7,5);
for i = [1,2,3,4,5,6]
    c1 = modelfunc1(beta3,i); c2 = modelfunc1(beta4,i);
    PRF_gx(i,:) = [0, c1, c2, c2, 2*pi];
    PRF_gy(i,:) = [0, -c1, -(c1/(2*pi - c1))*(2*pi - c2), 2*pi - c2, 0];
end
c1 = modelfunc1(beta3,100); c2 = modelfunc1(beta4,100);
PRF_gx(end,:) = [0, c1, c2, c2, 2*pi];
PRF_gy(end,:) = [0, -c1, -(c1/(2*pi - c1))*(2*pi - c2), 2*pi - c2, 0];

fig1h = figure(1);
set(fig1h, 'Position', [50, 200, 500, 340])
clf
hold on
plot(PRF_gx(1,:),PRF_gy(1,:), '-', 'Linewidth', 1.5, 'Color', [0,0,0])
plot(PRF_gx(2,:),PRF_gy(2,:), '-', 'Linewidth', 1.5, 'Color', [1,0,0])
plot(PRF_gx(3,:),PRF_gy(3,:), '-', 'Linewidth', 1.5, 'Color', [0,0.8,0])
plot(PRF_gx(4,:),PRF_gy(4,:), '-', 'Linewidth', 1.5, 'Color', [0,0,1])
plot(PRF_gx(5,:),PRF_gy(5,:), '-', 'Linewidth', 1.5, 'Color', [1,0.5,0.3])
plot(PRF_gx(6,:),PRF_gy(6,:), '-', 'Linewidth', 1.5, 'Color', [0.8,0.0,0.9])
plot(PRF_gx(7,:),PRF_gy(7,:), '-', 'Linewidth', 1.5, 'Color', [0.5,0.5,0.7])
hold off
grid on
%xlabelh = xlabel('Oscillator Phase, \phi');
%ylabelh = ylabel('General Phase Response, F_{g}(\phi)');
% For testing...
xlabelh = xlabel('Phase Value, \phi');
ylabelh = ylabel('Phase Response');

set(xlabelh,'Fontname','Times New Roman', 'Fontsize',12)
set(ylabelh,'Fontname','Times New Roman', 'Fontsize',12)
axis([0,2*pi,-pi-0.1,pi+0.1])
xticks(0:(pi/2):2*pi)
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
%yticks(-2*pi:(pi/2):2*pi)
%yticklabels({'-2\pi', '-3\pi/2', '-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
yticks([(-2*pi) (-7*pi/4) (-3*pi/2) (-5*pi/4) (-pi) (-3*pi/4) (-pi/2) (-pi/4) 0 (pi/4) (pi/2) (3*pi/4) pi (5*pi/4) (3*pi/2) (7*pi/4) 2*pi])
yticklabels({'-2\pi', '-7\pi/4', '-3\pi/2', '-5\pi/4', '-\pi', '-3\pi/4', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2', '7\pi/4', '2\pi'})

legendh = legend('$\delta^{-} = 1$', '$\delta^{-} = 2$', '$\delta^{-} = 3$', '$\delta^{-} = 4$',...
    '$\delta^{-} = 5$', '$\delta^{-} = 6$', '$\delta^{-} = \infty$', 'interpreter', 'latex');
set(legendh,'Fontname','Times New Roman', 'Fontsize',10, 'Location', 'northwest')
%}

%% Functions
function [lej_str] = LegendString(n,t)
if t == "all"
    t_new = "All-to-All";
elseif t == "ring"
    t_new = "Ring";
elseif t == "dring"
    t_new = "D-Ring";
elseif t == "lop3"
    t_new = "3 Ahead";
elseif t == "lop4"
    t_new = "4 Ahead";
elseif t == "near2far1"
    t_new = "Near2-Far1";
elseif t == "far3"
    t_new = "Farthest 3";
elseif t == "near4"
    t_new = "Nearest 4";
else
    t_new = "?";
end
lej_str = "N = " + string(n) + "; " + string(t_new); % Legend Info
end