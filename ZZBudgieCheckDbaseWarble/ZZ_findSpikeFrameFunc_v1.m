function spike_frame = ZZ_findSpikeFrameFunc_v1(spike_t, t)
%ZZ_FINDSPIKEFRAMEFUNC_V1 Summary of this function goes here
%   find what frames the spikes land
% spike_t is array of the spike time
% t is array of the spectrogram t
spike_frame = zeros(size(spike_t));
for spi=1:length(spike_t)
  for ti=2:length(t)
    if (t(ti-1)<=spike_t(spi)) && (t(ti)>=spike_t(spi))
      break
    end
  end
  spike_frame(spi) = ti;
end
end

