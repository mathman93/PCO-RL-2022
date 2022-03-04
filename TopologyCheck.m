% Check the topology of a given adjacency matrix
% Specifically, if it is connected or not

N = 11;
A_option = 'near4';
A = GetTopology(N,A_option);
%A = GetTopology(N,'custom',[2,4,5,7],true); % Custom topology

%% Specific ERG graph
N = 10;
A = zeros(N,N);
e1 = [1,1,1,2,2,2,3, 4,6,6,7,7, 8,8 ];
e2 = [2,5,6,3,4,8,10,8,7,8,9,10,9,10];
for k = 1:length(e1)
    i = e1(k); j = e2(k);
    A(i,j) = 1; A(j,i) = 1;
end

%%
CheckTopology(A);

xp = cos([0:N-1]*(2*pi/N));
yp = sin([0:N-1]*(2*pi/N));
nlabel10 = {'  1', '  2', '  3', '  4', '  5', '  6', '  7', '  8', '  9', '  10'};
nlabel11 = {'  1', '  2', '  3', '  4', '  5', '  6', '  7', '  8', '  9', '  10', '  11'};
% Figure 1: Plot of network topology
%
fig1 = figure(1);
clf
set(fig1, 'Position', [55, 50, 400, 400])
plt = plot(digraph(A),'o-', 'Layout', 'circle', 'NodeLabel', nlabel10);
set(plt, 'MarkerSize', 14, 'LineWidth', 1.5, 'NodeColor', [0.1,0.4,1], 'EdgeColor', [1,0.4,0], 'Linewidth', 2, 'EdgeAlpha', 0.5)
plt.XData = xp;
plt.YData = yp;
axis(1.15*[-1, 1, -1, 1])
%titleh = title("{``}Near4{''} Topology");
titleh = title("ERG Random Topology");
set(titleh,'Fontname','Times New Roman', 'Fontsize',13, 'interpreter', 'latex')
set(gca, 'box', 'off', 'xcolor', 'none', 'ycolor', 'none', 'Position', [0.03, 0.005, 0.94,0.92])
%}

%% Functions

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