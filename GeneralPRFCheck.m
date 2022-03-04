% Find optimal coefficients for general form of phase response function F_g
% from a given learned PRF

% Assume PRF is given (PRF1) from PCO_RL_v4.m (also PC_RL_v3.m)
% PRF1 is of size (P+1,N)

[c_matrix,c1,c2] = Optimize(PRF1,threshold);
disp([c1,c2]*(1/P))
disp([c1,c2]*(threshold/P))


%% Graphs
% Best fit measures
%{
figure(1)
clf
imagesc(c_matrix)
colorbar
caxis([min(min(c_matrix)), min(min(c_matrix))+100])
%}

%% Functions
function [c_matrix, c1_index, c2_index] = Optimize(PRF1,threshold)
[P_prime, N] = size(PRF1);
P = P_prime - 1;
c1_index = zeros(1,N); % final learned parameter values
c2_index = zeros(1,N);
% Perform first approximation of parameter values
c_matrix = 1000*ones(P+1,P+1); % Initial values
for n = 1:N % For each oscillator
    % Check each possible comination of c1 and c2
    for c2 = 0:P
        for c1 = 0:c2 % Ensure c1 <= c2
            % Compute differences between PRF1 and GeneralPRF
            diffs = zeros(1,P+1);
            for i = 1:(P+1) % For each data point in PRF1
                if i<c2 && i>=c1
                    scale = 1.2; % Scale factor for center portion of GeneralPRF
                else
                    scale = 1.0; % Scale factor for end portions of GeneralPRF
                end
                diffs(i) = scale * (PRF1(i,n) - GeneralPRF((i-1)*(threshold/P), threshold, c1*(1/P), c2*(1/P)));
            end
            %diff_sum = sum(diffs);
            c_matrix(c1+1,c2+1) = sum(diffs.^2); % Error = Sum of square of differences
            %c_matrix(c1+1,c2+1) = sum(abs(diffs));
        end
    end
    [c2_min, c2_i] = min(c_matrix,[],2);
    [~, c1_i] = min(c2_min);
    %min(min(c_matrix)); % Minimum value
    %c_matrix(c1_i,c2_i(c1_i)); % Alternate expression of minimum value
    c1_index(n) = c1_i-1;
    c2_index(n) = c2_i(c1_i)-1;
end
end

function [c1_value, c2_value] = Optimize2(PRF1, threshold, width)
[P_prime, N] = size(PRF1);
P = P_prime - 1;
nu = 1/width;
K = 2*nu + 1;
c1_cent = c1_index*(1/P);
c2_cent = c2_index*(1/P);
end

% Third PRF for PCO Synchronization (based on RL training)
% Parameters:
%   phase = float; single oscillator phase
%   threshold = float; oscillator threshold
%   refractory = float; PRF refractory period
% Returns:
%   y = float; the value of the PRF at phase
function [y] = GeneralPRF(phase, threshold, c1, c2)
%c1 = 0.35; % [0.3, 0.45]
%c2 = 0.75; % [0.7, 0.8]
if phase <= threshold
    if phase < c1 * threshold
        y = -phase;
    elseif phase < c2 * threshold
        y = -(c1/(1-c1))*(threshold - phase);
    else % c2 * threshold < phase < threshold
        y = threshold - phase;
    end
else
    y = 0.0;
end
end