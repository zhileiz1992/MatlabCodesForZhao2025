function [diffAmp] = ZZ_CalculateFreqDiff(signal, fs, ampIgnore, freqRange1, freqRange2)
% calculate a difference in amplitude between freqRange1 and freqRange2
addpath(genpath("\\goldbergsongserver2.nbb.cornell.edu\z4\zz367\BudgieContactCalls\FromHan"));
addpath(genpath("\\goldbergsongserver2.nbb.cornell.edu\z4\zz367\BudgieContactCalls\BudgieContactCallsAnalysis")); 

%% band pass signal
highPass = 500;
lowPass = 10000; 
[b,a] = butter(2, highPass/fs, 'high');
signalFiltered = filtfilt(b, a, signal);
[b,a] = butter(2, lowPass/fs, 'low');
signalFiltered = filtfilt(b, a, signalFiltered);

%% calculate spectrogram
[spec, dt, f ,T] = get_spec(signalFiltered, fs);
% ignore where amplitude is very small
spec(spec<ampIgnore) = ampIgnore;

%% calculate amp difference between two freq ranges
idxLowFreq = (f>=freqRange1(1)) & (f<=freqRange1(2));
idxOtherFreq = (f>=freqRange2(1)) & (f<=freqRange2(2));
% calculate an average power
lowFreqAmp = sum(spec(idxLowFreq,:), 1)/sum(idxLowFreq);
otherFreqAmp = sum(spec(idxOtherFreq,:), 1)/sum(idxOtherFreq);
% smoothing
lowFreqAmp = smooth(lowFreqAmp, 10);
otherFreqAmp = smooth(otherFreqAmp, 10);
% figure; plot(lowFreqAmp, 'b');
% hold on; plot(otherFreqAmp, 'g');
% calculat the difference
diffAmp = otherFreqAmp - lowFreqAmp;
% hold on; plot(diffAmp, 'r');
end

