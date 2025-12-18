% example script from GPT
% Generate example data
rng(1); % for reproducibility
X = [randn(100,2)*0.75 + ones(100,2);
     randn(100,2)*0.5 - ones(100,2);
     randn(100,2)*0.6 + [2 -2]];

% Choose number of clusters
k = 3;

% Perform k-means clustering
[idx, C] = kmeans(X, k, 'Replicates', 10);

% Silhouette analysis
figure;
silhouette(X, idx);
title(sprintf('Silhouette Plot for k = %d', k));
xlabel('Silhouette Value');
ylabel('Cluster');


k_list = 2:6;
meanSilh = zeros(size(k_list));

for i = 1:length(k_list)
    k = k_list(i);
    idx = kmeans(X, k, 'Replicates', 10);
    s = silhouette(X, idx);
    meanSilh(i) = mean(s);
end

figure;
plot(k_list, meanSilh, '-o');
xlabel('Number of Clusters (k)');
ylabel('Mean Silhouette Value');
title('Silhouette Analysis for Different k');