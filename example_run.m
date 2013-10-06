% Demonstration of the use of RW_ModMax.
% Based on "maindriver.m" by Lei Tang (see http://leitang.net/code/social-dimension/SocioDim.zip)
% Robin Devooght 2013, october 4th

global network;
load blogcatalog.mat



% Extract top eigenvectors of the random walk based modularity matrix for use as nodes' features
options.k = 500; % Number of eigenvectors
options.theta = 1; % Inverse temperature
eigenvectors = RW_ModMax(options); 


% randomly generate index_tr (training nodes index) and index_te (test nodes index)
n = size(network, 1);
index = randperm(n);
index_tr = index(1:ceil(0.1*n));  % 10% labeled nodes for training
index_te = index(1+ceil(0.1*n):end);  % 90% unlabeled nodes for test
labels = group(index_tr, :); % the labels of nodes for training

% build the classifier and make predictions
C = 20; % the C parameter in SVM Classifier
[predscore] = SocioDim(eigenvectors, labels, index_tr, index_te, C);

[perf, pred] = evaluate(predscore, group(index_te, :));

perf.micro_F1
perf.macro_F1
perf.acc