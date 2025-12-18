function [dbase, dbaseOld] = ZZ_flushDbaseSegmentationFunc_v1(dbase, onsets, offsets)
% flush syllable segmentation into a dbase
% assume that each cell in onsets/offsets correspond to one audio file in dbase
dbaseOld = dbase;

% only change the segmentation part of dbase
% leave other things untouched
for fi = 1:length(dbase.SoundFiles)
%   fi = 51;
  seg_times = [onsets{fi}; offsets{fi}];
  seg_titles = {};
  seg_selected = [];
  for ii=1:length(onsets{fi})
    seg_titles{ii} = 'a';
    seg_selected(ii) = 0;
  end
  dbase.SegmentThresholds(fi) = 1;
  dbase.SegmentTimes{fi} = seg_times';
  dbase.SegmentTitles{fi} = seg_titles;
  dbase.SegmentIsSelected{fi} = seg_selected;
end

end

