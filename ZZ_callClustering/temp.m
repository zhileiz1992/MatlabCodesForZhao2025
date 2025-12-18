function d = emd_1d(p, q)
    d = sum(abs(cumsum(p) - cumsum(q)));
end
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