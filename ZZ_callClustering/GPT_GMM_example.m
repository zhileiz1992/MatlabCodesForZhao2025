%% Synthetic data
rng(1); % For reproducibility
X_train = [mvnrnd([2 2], [1 0.8; 0.8 1], 100);
           mvnrnd([-2 -3], [1 -0.5; -0.5 1], 100);
           mvnrnd([5 -2], [1 0; 0 0.5], 100)];

% Specify number of components (clusters)
K = 3;

% Fit Gaussian Mixture Model with full covariance
options = statset('MaxIter', 500); % increase iterations if needed
gmm = fitgmdist(X_train, K, ...
    'CovarianceType', 'full', ...
    'RegularizationValue', 1e-6, ... % for numerical stability
    'Options', options, ...
    'Replicates', 10); % run EM multiple times to avoid local optima

% Predict cluster labels
idx = cluster(gmm, X_train); % assign each point to a component

% Plot clustering result
figure;
gscatter(X_train(:,1), X_train(:,2), idx);
hold on;
% Plot GMM component means
plot(gmm.mu(:,1), gmm.mu(:,2), 'kx', 'LineWidth', 2, 'MarkerSize', 10);
title('GMM Clustering with Full Covariance');
xlabel('Feature 1'); ylabel('Feature 2');
legend('Cluster 1','Cluster 2','Cluster 3','GMM Means');


%% Real Iris dataset
% Load the Iris dataset
clear; close all;
load fisheriris % variables: meas (150x4), species (cell array)
X = meas; % Use the 4-dimensional feature vectors

% Number of clusters (we know there are 3 Iris species)
K = 3;

% Fit a GMM with full covariance using Expectation-Maximization
options = statset('MaxIter', 500);
gmm = fitgmdist(X, K, ...
    'CovarianceType', 'full', ...
    'RegularizationValue', 1e-6, ...
    'Replicates', 10, ...
    'Options', options);

% Assign clusters to each observation
% idx = cluster(gmm, X);
[idx,NLOGL,POST] = cluster(gmm, X);

% Compare GMM clustering with actual species
figure;
gscatter(X(:,1), X(:,2), idx);
xlabel('Sepal Length');
ylabel('Sepal Width');
title('GMM Clustering of Iris Dataset (First 2 Features)');
legend('Cluster 1','Cluster 2','Cluster 3','Location','best');

% Optional: silhouette analysis
figure;
silhouette(X, idx);
title('Silhouette Plot for GMM Clustering (Iris Data)');

