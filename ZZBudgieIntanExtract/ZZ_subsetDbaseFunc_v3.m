function [dbase_subset, dbase_old] = ZZ_subsetDbaseFunc_v3(dbase, select_idx)
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

% to work with Brian's new gui
dbase_subset.Notes = {};
for ii=1:length(select_idx)
  dbase_subset.Notes{ii} = dbase_old.Notes{select_idx(ii)};
end
temp = dbase_old.FileReadState;
dbase_subset.FileReadState = temp(select_idx);

if isfield(dbase, 'Properties')
  a = dbase_old.Properties;
  if isempty(a)
    dbase_subset.Properties = a;
  else
    dbase_subset.Properties = a(select_idx,:);
  end
end

end

