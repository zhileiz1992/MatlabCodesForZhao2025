function [dbase_subset, dbase_old] = ZZ_subsetDbaseFunc_v1(dbase, select_idx)
% subset a dbase
% create a new base that only have the select_idx in it
dbase_old = dbase;
dbase_subset = dbase;
rep_fields = {'Times', 'FileLength', 'SoundFiles', 'SegmentThresholds', 'SegmentTimes', 'SegmentTitles', 'SegmentIsSelected', ...
              'MarkerTimes', 'MarkerTitles', 'MarkerIsSelected'};
for ri=1:length(rep_fields)
  rf = rep_fields{ri};
  a = dbase.(rf);
  dbase_subset.(rf) = a(select_idx);
end
% deal with the channel files
for ch=1:length(dbase.ChannelFiles)
  a = dbase.ChannelFiles{ch};
  dbase_subset.ChannelFiles{ch} = a(select_idx);
end

end

