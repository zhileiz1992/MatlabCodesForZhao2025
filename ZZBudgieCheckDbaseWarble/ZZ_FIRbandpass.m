function [signalFiltered] = ZZ_FIRbandpass(signal, fs, freq1, freq2, ord)
  b = fir1(ord,[freq1 freq2]/(fs/2));
  signalFiltered = filtfilt(b, 1, signal);
end

