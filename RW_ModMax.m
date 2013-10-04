function eigenvectors = RW_ModMax(options)
%
% Compute main eigenvectors of the random walk-based modularity matrix without computing explicitly the modularity matrix itself.
% Adapted for large scale networks.
%
% INPUT :
% - options
%   - options.theta > 0 tune the respective influence of short and long
%     paths in the bag-of-paths model. default: 1.
%   - options.features : number of features to extract from BoP modularity
%     matrix. default: 100.
% NB: the network adjacency matrix is the global variable "network".
%
% OUTPUT : 
% - eigenvectors is a n*k matrix where n is the number of nodes in the graph
%   and k is given in input. It contains the k dominant eigenvectors of the
%   BoP modularity matrix.
%
% Robin Devooght : 2013, october 4th

global network;

% Gives default values to unspecified options
if ~isfield(options, 'theta')
    options.theta = 1;
end
if ~isfield(options, 'features')
    options.features = 100;
end

[n, m] = size(network);

if n ~= m
    error('RW_ModMax:A_square', 'Adjacency matrix must be square');
end
if options.theta <= 0
    error('RW_ModMax:theta_range', 'theta must be positive');
end
if options.features <= 0 || options.features > n
    error('RW_ModMax:k_range', 'the number of features must be comprised between 1 and n (the number of nodes in the graph)');
end

%% Compute the W matrix
nonZeroIndex = network > 0;

Cost = sparse(n,n);
Cost(nonZeroIndex) = 1./(network(nonZeroIndex));
D = diag(sum(network,2).^(-1));
Pref = D*network;

W = sparse(n,n);
W(nonZeroIndex) = exp(-options.theta*Cost(nonZeroIndex)).*Pref(nonZeroIndex);

%% Computation of various quantities needed in the "modularity_product" function.
e = ones(n,1);
Ze = Zx(W, e, e); % estimate Z*e, with Z = inv(I - W)
z = sum(Ze); % equal to sum(sum(Z)) where Z is the fondamental matrix : Z = inv(I-W)

%% Extract dominant eigen vectors using "eigs" routine. 
%"eigs" relies on "modularity_product" to compute the product of the modularity matrix with a column vector.
opts.issym = 1;
opts.isreal = 1;

[eigenvectors, D] = eigs(@(x)modularity_product(x), n, options.features, 'LA', opts);


function [Q_times_x] = modularity_product(x)
    % Compute Q*x where Q is the modularity matrix, without computing Q.
    y = Zx(W, x, e);

    Q_times_x = y/z - sum(y)*Ze/z^2;
end
end
