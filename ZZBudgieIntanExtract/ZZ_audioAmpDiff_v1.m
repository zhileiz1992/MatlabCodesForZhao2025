% function to calculate amplitude diffeerence between two audio streams
function [amp_diff, delta_t_amp, pos_idx, neg_idx] = ZZ_audioAmpDiff_v1(signal1, signal2, fs, plot_bool)
% band-filter data
highPass = 500;
lowPass = 10000;
signal1 = ZZ_bandPass_v1(signal1', fs, highPass, lowPass);
signal2 = ZZ_bandPass_v1(signal2', fs, highPass, lowPass);

% calculate amplitude
signalComb = [signal1 signal2];
amp = acousticLoudness(signalComb,fs,'SoundField','diffuse','TimeVarying',true, 'TimeResolution', 'standard');
% smooth 
amp_smooth1 = smooth(amp(:,1), 50);
amp_smooth2 = smooth(amp(:,2), 50);
% calculate diff 
amp_diff = amp_smooth1 - amp_smooth2;
% standard resolution is 2 ms
delta_t_amp = 0.002; 

% also infer who is the singer
% based on the number of points that amp_diff cross a threshold
thre = 0.5; 
pos_idx = find(amp_diff>=thre);
neg_idx = find(amp_diff<=-thre);

if plot_bool
  figure; 
  subplot(3,1,1); plot(amp);
  subplot(3,1,2); plot(amp_smooth1); hold on; plot(amp_smooth2);
  subplot(3,1,3); plot(amp_diff); hold on; yline(0);
  hold on; plot(pos_idx, amp_diff(pos_idx), 'r.');
  hold on; plot(neg_idx, amp_diff(neg_idx), 'g.');
end

end