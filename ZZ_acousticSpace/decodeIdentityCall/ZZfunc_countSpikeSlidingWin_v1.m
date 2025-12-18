function output = ZZfunc_countSpikeSlidingWin_v1(aligned_spike, bin_size, hop_size)
    % Get dimensions
    [n, m] = size(aligned_spike);
    
    % Calculate the number of sliding windows
    num_windows = floor((m - bin_size) / hop_size) + 1;
    
    % Initialize output array
    output = zeros(n, num_windows);
    
    % Loop over each trial
    for i = 1:n
        % Loop over each window
        for j = 1:num_windows
            start_idx = 1 + (j - 1) * hop_size;
            end_idx = start_idx + bin_size - 1;
            output(i, j) = sum(aligned_spike(i, start_idx:end_idx));
        end
    end
end
