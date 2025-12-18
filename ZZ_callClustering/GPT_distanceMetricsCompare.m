% a script to compare different distance metrics
clear; close all

%% same freq-time pattern, different amplitude
% Parameters
fs = 8000;                  % Sampling frequency
t = 0:1/fs:1-1/fs;          % 1 second duration

% Generate two signals with same spectral content but different loudness
signalA = sin(2*pi*440*t) + 0.5*sin(2*pi*880*t);
signalB = 3 * signalA;      % Louder but same shape

% Compute spectrograms
win = hamming(256);
[S_A, F, T] = spectrogram(signalA, win, 128, 256, fs);
[S_B, ~, ~] = spectrogram(signalB, win, 128, 256, fs);

% Log-magnitude spectrograms
log_SA = log1p(abs(S_A));
log_SB = log1p(abs(S_B));

% Flatten for distance calculations
vecA = reshape(log_SA, [], 1);
vecB = reshape(log_SB, [], 1);

% Normalize and center
vecA_norm = vecA / norm(vecA);
vecB_norm = vecB / norm(vecB);
vecA_centered = vecA - mean(vecA);
vecB_centered = vecB - mean(vecB);

% Compute distances
cosine_dist = 1 - dot(vecA_norm, vecB_norm);
correlation_dist = 1 - dot(vecA_centered, vecB_centered) / ...
                   (norm(vecA_centered) * norm(vecB_centered));
euclidean_dist = norm(vecA - vecB);

fprintf('Cosine distance: %.6f\n', cosine_dist);
fprintf('Correlation distance: %.6f\n', correlation_dist);
fprintf('Euclidean distance: %.6f\n', euclidean_dist);

% Set common color limits for consistent visualization
clim = [min([log_SA(:); log_SB(:)]), max([log_SA(:); log_SB(:)])];

% Plot both spectrograms
figure;

subplot(1,2,1);
imagesc(T, F, log_SA, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Log-Spectrogram: Signal A');
colorbar;

subplot(1,2,2);
imagesc(T, F, log_SB, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Log-Spectrogram: Signal B');
colorbar;


%% same amplitude, frequecy slight shifted 100hz
% Parameters
fs = 8000;                  % Sampling frequency
t = 0:1/fs:1-1/fs;          % 1 second duration

% Generate signals
signalA = sin(2*pi*440*t) + 0.5*sin(2*pi*880*t);        % base signal
signalB = sin(2*pi*540*t) + 0.5*sin(2*pi*980*t);        % frequency-shifted by +100 Hz

% Spectrogram computation
win = hamming(256);
[S_A, F, T] = spectrogram(signalA, win, 128, 256, fs);
[S_B, ~, ~] = spectrogram(signalB, win, 128, 256, fs);

% Log-magnitude
log_SA = log1p(abs(S_A));
log_SB = log1p(abs(S_B));

% Flatten into vectors
vecA = reshape(log_SA, [], 1);
vecB = reshape(log_SB, [], 1);

% Normalize and center
vecA_norm = vecA / norm(vecA);
vecB_norm = vecB / norm(vecB);
vecA_centered = vecA - mean(vecA);
vecB_centered = vecB - mean(vecB);

% Compute distances
cosine_dist = 1 - dot(vecA_norm, vecB_norm);
correlation_dist = 1 - dot(vecA_centered, vecB_centered) / ...
                   (norm(vecA_centered) * norm(vecB_centered));
euclidean_dist = norm(vecA - vecB);

fprintf('Cosine distance: %.6f\n', cosine_dist);
fprintf('Correlation distance: %.6f\n', correlation_dist);
fprintf('Euclidean distance: %.6f\n', euclidean_dist);

% Common color limits
clim = [min([log_SA(:); log_SB(:)]), max([log_SA(:); log_SB(:)])];

% Plot spectrograms
figure;

subplot(1,2,1);
imagesc(T, F, log_SA, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Log-Spectrogram: Signal A');
colorbar;

subplot(1,2,2);
imagesc(T, F, log_SB, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Log-Spectrogram: Signal B (+100 Hz)');
colorbar;


%% Shirt the frequency band by 100ms
% Parameters
fs = 8000;                  % Sampling frequency
t = 0:1/fs:1-1/fs;          % 1 second duration

% Create a frequency sweep (chirp) that peaks in the middle
f0 = 400;                   % starting frequency
f1 = 1000;                  % peak frequency
t_peak = 0.5;               % peak at 0.5s

% Signal A: a Gaussian-modulated chirp peaking at 0.5s
sweepA = chirp(t, f0, t_peak, f1, 'linear');
envA = exp(-((t - t_peak).^2) / (2*(0.05)^2));  % narrow Gaussian envelope
signalA = sweepA .* envA;

% Signal B: same sweep but delayed by 100 ms
delay = 0.1;  % 100 ms
tB = t - delay;
sweepB = chirp(tB, f0, t_peak, f1, 'linear');
envB = exp(-((tB - t_peak).^2) / (2*(0.05)^2));
signalB = sweepB .* envB;
signalB(tB < 0) = 0;  % zero-pad beginning

% Compute spectrograms
win = hamming(256);
[S_A, F, T] = spectrogram(signalA, win, 128, 256, fs);
[S_B, ~, ~] = spectrogram(signalB, win, 128, 256, fs);

% Log-magnitude
log_SA = log1p(abs(S_A));
log_SB = log1p(abs(S_B));

% Flatten to 1D
vecA = reshape(log_SA, [], 1);
vecB = reshape(log_SB, [], 1);

% Normalize and center
vecA_norm = vecA / norm(vecA);
vecB_norm = vecB / norm(vecB);
vecA_centered = vecA - mean(vecA);
vecB_centered = vecB - mean(vecB);

% Compute distances
cosine_dist = 1 - dot(vecA_norm, vecB_norm);
correlation_dist = 1 - dot(vecA_centered, vecB_centered) / ...
                   (norm(vecA_centered) * norm(vecB_centered));
euclidean_dist = norm(vecA - vecB);

fprintf('Cosine distance: %.6f\n', cosine_dist);
fprintf('Correlation distance: %.6f\n', correlation_dist);
fprintf('Euclidean distance: %.6f\n', euclidean_dist);

% Common color limits
clim = [min([log_SA(:); log_SB(:)]), max([log_SA(:); log_SB(:)])];

% Plot
figure;

subplot(1,2,1);
imagesc(T, F, log_SA, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Signal A (Peak @ 0.5s)');
colorbar;

subplot(1,2,2);
imagesc(T, F, log_SB, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Signal B (Peak @ 0.6s)');
colorbar;


%% Different freq-time pattern
% Parameters
fs = 8000;
t = 0:1/fs:1-1/fs;

% Envelope (same for both)
t_peak = 0.5;
env = exp(-((t - t_peak).^2) / (2*(0.05)^2));  % Gaussian envelope

% Signal A: upward chirp
sweepA = chirp(t, 400, t_peak, 1000, 'linear');
signalA = sweepA .* env;

% Signal B: downward chirp
sweepB = chirp(t, 1000, t_peak, 400, 'linear');
signalB = sweepB .* env;

% Spectrograms
win = hamming(256);
[S_A, F, T] = spectrogram(signalA, win, 128, 256, fs);
[S_B, ~, ~] = spectrogram(signalB, win, 128, 256, fs);

% Log-magnitude
log_SA = log1p(abs(S_A));
log_SB = log1p(abs(S_B));

% Flatten to 1D vectors
vecA = reshape(log_SA, [], 1);
vecB = reshape(log_SB, [], 1);

% Normalize and center
vecA_norm = vecA / norm(vecA);
vecB_norm = vecB / norm(vecB);
vecA_centered = vecA - mean(vecA);
vecB_centered = vecB - mean(vecB);

% Compute distances
cosine_dist = 1 - dot(vecA_norm, vecB_norm);
correlation_dist = 1 - dot(vecA_centered, vecB_centered) / ...
                   (norm(vecA_centered) * norm(vecB_centered));
euclidean_dist = norm(vecA - vecB);

fprintf('Cosine distance: %.6f\n', cosine_dist);
fprintf('Correlation distance: %.6f\n', correlation_dist);
fprintf('Euclidean distance: %.6f\n', euclidean_dist);

% Set common color scale
clim = [min([log_SA(:); log_SB(:)]), max([log_SA(:); log_SB(:)])];

% Plot
figure;

subplot(1,2,1);
imagesc(T, F, log_SA, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Log-Spectrogram: Upward Chirp (A)');
colorbar;

subplot(1,2,2);
imagesc(T, F, log_SB, clim); axis xy;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Log-Spectrogram: Downward Chirp (B)');
colorbar;


%% compare raw with time wrapped
% Parameters
fs = 8000;
t = 0:1/fs:1-1/fs;
t_peak = 0.5;
env = exp(-((t - t_peak).^2) / (2*(0.05)^2));  % Gaussian envelope

% Base signal: upward chirp
signal_base = chirp(t, 400, t_peak, 1000, 'linear') .* env;

% Time-shifted version (e.g., 20 ms delay)
delay_samples = round(0.020 * fs);  % 20 ms
signal_time_shift = [zeros(1, delay_samples), signal_base(1:end-delay_samples)];

% Frequency-shifted version (e.g., +50 Hz)
signal_freq_shift = chirp(t, 450, t_peak, 1050, 'linear') .* env;

% Spectrograms
win = hamming(256);
[S_base, F, T] = spectrogram(signal_base, win, 128, 256, fs);
[S_time, ~, ~] = spectrogram(signal_time_shift, win, 128, 256, fs);
[S_freq, ~, ~] = spectrogram(signal_freq_shift, win, 128, 256, fs);

% Log-magnitude
log_base = log1p(abs(S_base));
log_time = log1p(abs(S_time));
log_freq = log1p(abs(S_freq));

% Vectorize
vec_base = reshape(log_base, [], 1);
vec_time = reshape(log_time, [], 1);
vec_freq = reshape(log_freq, [], 1);

% Distance function
compute_dists = @(a, b) struct( ...
    'cosine', 1 - dot(a/norm(a), b/norm(b)), ...
    'correlation', 1 - corr(a, b), ...
    'euclidean', norm(a - b));

% Compute distances
dist_time = compute_dists(vec_base, vec_time);
dist_freq = compute_dists(vec_base, vec_freq);

% Display
fprintf('Distance to Time-Shifted Signal (20 ms):\n');
fprintf('  Cosine:      %.6f\n', dist_time.cosine);
fprintf('  Correlation: %.6f\n', dist_time.correlation);
fprintf('  Euclidean:   %.6f\n\n', dist_time.euclidean);

fprintf('Distance to Frequency-Shifted Signal (+50 Hz):\n');
fprintf('  Cosine:      %.6f\n', dist_freq.cosine);
fprintf('  Correlation: %.6f\n', dist_freq.correlation);
fprintf('  Euclidean:   %.6f\n', dist_freq.euclidean);

% Common color range
clim = [min([log_base(:); log_time(:); log_freq(:)]), ...
        max([log_base(:); log_time(:); log_freq(:)])];

% Plot
figure;
subplot(1,3,1);
imagesc(T, F, log_base, clim); axis xy;
xlabel('Time (s)'); ylabel('Freq (Hz)');
title('Base Signal');
colorbar;

subplot(1,3,2);
imagesc(T, F, log_time, clim); axis xy;
xlabel('Time (s)'); ylabel('Freq (Hz)');
title('Time-Shifted (20 ms)');
colorbar;

subplot(1,3,3);
imagesc(T, F, log_freq, clim); axis xy;
xlabel('Time (s)'); ylabel('Freq (Hz)');
title('Freq-Shifted (+50 Hz)');
colorbar;


%% wrapped version
% Signal setup (same as before)
% Signal setup (same as before)
fs = 8000;
t = 0:1/fs:1-1/fs;
t_peak = 0.5;
env = exp(-((t - t_peak).^2) / (2*(0.05)^2));

signal_base = chirp(t, 400, t_peak, 1000, 'linear') .* env;
delay_samples = round(0.020 * fs);
signal_time_shift = [zeros(1, delay_samples), signal_base(1:end-delay_samples)];
signal_freq_shift = chirp(t, 450, t_peak, 1050, 'linear') .* env;

% Spectrograms
win = hamming(256);
[S_base, F, T] = spectrogram(signal_base, win, 128, 256, fs);
[S_time, ~, ~] = spectrogram(signal_time_shift, win, 128, 256, fs);
[S_freq, ~, ~] = spectrogram(signal_freq_shift, win, 128, 256, fs);

log_base = log1p(abs(S_base));
log_time = log1p(abs(S_time));
log_freq = log1p(abs(S_freq));

% ---------- 1. DTW for Time Shift ----------
% Use mean frequency across time bins as time series
mean_base = mean(log_base, 1);   % [1 x T]
mean_time = mean(log_time, 1);   % [1 x T]

% DTW alignment (requires Signal Processing Toolbox)
[~, ix_time, ix_base] = dtw(mean_time, mean_base);

% Warp log_time to match log_base's time axis
log_time_dtw = log_time(:, ix_time);  % aligned version

% Pad to same length
min_len = min(size(log_base,2), size(log_time_dtw,2));
vec_base = reshape(log_base(:,1:min_len), [], 1);
vec_time_dtw = reshape(log_time_dtw(:,1:min_len), [], 1);

% ---------- 2. Spectral Cross-Correlation for Frequency Shift ----------
% Cross-correlate average spectra to find best frequency shift
mean_spec_base = mean(log_base, 2);  % [F x 1]
mean_spec_freq = mean(log_freq, 2);  % [F x 1]
[cross_corr, lags] = xcorr(mean_spec_freq, mean_spec_base, 'coeff');
[~, max_idx] = max(cross_corr);
lag_best = lags(max_idx);  % shift in frequency bins

% Circularly shift the spectrogram to align
log_freq_shifted = circshift(log_freq, [-lag_best, 0]);

% Truncate to match base size
vec_freq_shifted = reshape(log_freq_shifted, [], 1);
vec_base_full = reshape(log_base, [], 1);

% ---------- Compute Distances ----------
compute_dists = @(a, b) struct( ...
    'cosine', 1 - dot(a/norm(a), b/norm(b)), ...
    'correlation', 1 - corr(a, b), ...
    'euclidean', norm(a - b));

% Before vs after DTW (time alignment)
vec_time_unaligned = reshape(log_time(:,1:min_len), [], 1);
d_time_before = compute_dists(vec_base, vec_time_unaligned);
d_time_after  = compute_dists(vec_base, vec_time_dtw);

% Before vs after spectral alignment
d_freq_before = compute_dists(vec_base_full, reshape(log_freq, [], 1));
d_freq_after  = compute_dists(vec_base_full, vec_freq_shifted);

% ---------- Print Results ----------
fprintf('\n=== Time-Shift (20 ms) ===\n');
fprintf('Before DTW:  Cosine: %.4f  Corr: %.4f  Euc: %.4f\n', ...
    d_time_before.cosine, d_time_before.correlation, d_time_before.euclidean);
fprintf('After  DTW:  Cosine: %.4f  Corr: %.4f  Euc: %.4f\n', ...
    d_time_after.cosine, d_time_after.correlation, d_time_after.euclidean);

fprintf('\n=== Frequency-Shift (+50 Hz) ===\n');
fprintf('Before Align: Cosine: %.4f  Corr: %.4f  Euc: %.4f\n', ...
    d_freq_before.cosine, d_freq_before.correlation, d_freq_before.euclidean);
fprintf('After  Align: Cosine: %.4f  Corr: %.4f  Euc: %.4f\n', ...
    d_freq_after.cosine, d_freq_after.correlation, d_freq_after.euclidean);
  
  
%% Earth Mover Distance (EMD)
% Parameters
fs = 8000;
t = 0:1/fs:1-1/fs;
t_peak = 0.5;
env = exp(-((t - t_peak).^2) / (2*(0.05)^2));

% Signal A: chirp from 400 to 1000 Hz
signalA = chirp(t, 400, t_peak, 1000, 'linear') .* env;

% Signal B: chirp from 500 to 1100 Hz (100 Hz upward shift)
signalB = chirp(t, 500, t_peak, 1100, 'linear') .* env;

% Compute spectrograms
win = hamming(256);
[S_A, F, ~] = spectrogram(signalA, win, 128, 256, fs);
[S_B, ~, ~] = spectrogram(signalB, win, 128, 256, fs);

% Take mean power across time (1D histogram of frequency energy)
histA = mean(abs(S_A), 2);
histB = mean(abs(S_B), 2);

% Normalize to probability distributions
pA = histA / sum(histA);
pB = histB / sum(histB);

% Define ground distance: pairwise frequency bin distances (Euclidean)
D = pdist2(F, F);  % F is frequency vector

% Compute EMD (Computer Vision Toolbox required)
% emd_val = emd(pA, pB, D);
emd_val = emd_1d(pA, pB);  

% Also compute cosine and Euclidean distances
cosine_val = 1 - dot(pA, pB) / (norm(pA) * norm(pB));
euclidean_val = norm(pA - pB);

% Print results
fprintf('=== Distance Between Mean Spectra ===\n');
fprintf('EMD:       %.6f\n', emd_val);
fprintf('Cosine:    %.6f\n', cosine_val);
fprintf('Euclidean: %.6f\n', euclidean_val);

% Plot spectra
figure;
plot(F, pA, 'b-', 'LineWidth', 1.5); hold on;
plot(F, pB, 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
legend('Signal A', 'Signal B');
title('Mean Spectra Comparison');
grid on;

