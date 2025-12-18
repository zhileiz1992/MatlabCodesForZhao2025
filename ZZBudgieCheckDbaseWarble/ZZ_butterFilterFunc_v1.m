function [signalFiltered] = ZZ_butterFilterFunc_v1(signal, fs, highPass, lowPass, order)
%ZZ_BUTTERFILTERFUNC_V1 Summary of this function goes here
% band-pass filtered raw signal
[b,a] = butter(order, highPass/fs, 'high');
signalFiltered = filtfilt(b, a, signal);
[b,a] = butter(order, lowPass/fs, 'low');
signalFiltered = filtfilt(b, a, signalFiltered);
end

