function [signalFiltered] = ZZ_bandPass_v1(signal, fs, highPass, lowPass)
%% band pass signal
[b,a] = butter(2, highPass/fs, 'high');
signalFiltered = filtfilt(b, a, signal);
[b,a] = butter(2, lowPass/fs, 'low');
signalFiltered = filtfilt(b, a, signalFiltered);

end

