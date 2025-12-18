function [flatness, dt] = ZZ_CalculateFlatness(signal, fs, ampIgnore, maskFreq)
% calculate flatness of sound signal, adapted from Han's script

%% band pass signal
highPass = 500;
lowPass = 10000; 
[b,a] = butter(2, highPass/fs, 'high');
signalFiltered = filtfilt(b, a, signal);
[b,a] = butter(2, lowPass/fs, 'low');
signalFiltered = filtfilt(b, a, signalFiltered);

%% calculate tonal flatness
[spec, dt, f ,T] = get_spec(signalFiltered, fs);
% ignore where amplitude is very small
spec(spec<ampIgnore) = ampIgnore;
flatness=(geomean(exp(spec(f>maskFreq,:)),1)./nanmean(exp(spec(f>maskFreq,:)),1));
% smoothing 
flatness = smooth(flatness,10);
end

