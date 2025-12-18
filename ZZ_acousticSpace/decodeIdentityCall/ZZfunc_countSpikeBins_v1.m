function [bin_sums] = ZZfunc_countSpikeBins_v1(aligned_spike, bin_pt)
% count the number of spikes in non-overlapping bins
[n, m] = size(aligned_spike); % Get dimensions of input array
num_bins = floor(m / bin_pt); % Number of complete bins per row

% Initialize output array: n rows, num_bins columns
bin_sums = zeros(n, num_bins);

% Process each row
for i = 1:n
  % Truncate row to fit complete bins
  row_data = aligned_spike(i, 1:num_bins * bin_pt);
  % Reshape into a matrix where each column is a bin
  binned_data = reshape(row_data, bin_pt, num_bins);
  % Sum each bin (column) and store in output
  bin_sums(i, :) = sum(binned_data, 1);
end

end

