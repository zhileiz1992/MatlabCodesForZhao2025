figure; 
plot(flatness);
% add the threshold lines
yline(-param.thresholdFlatness, 'g'); hold on;
ylim([-0.05 1.05]);
% mark the onsets and offsets of identified syllables


fd = '/mnt/z4/zz367/EphysMONAO/Analyzed/ExtractedIntan/pair1/2023-11-04/pair1CU21RigB/NCFiles/warble';
for ch=0:21
  fn = fullfile(fd, sprintf('warble_00150_00006_1_pair1CU21_231104_230211_chan%d.nc', ch));
  data = ncread(fn, 'data');
  disp(size(data));
end

% check if ephys files are correctly parsed
ri=150;
% read data in
dataAll = [];
fns_rhd = metaInfo(ri).filenames;
for fi=1:length(fns_rhd)
  data = read_Intan_RHD2000_file_to_struct_2_ZZ(fns_rhd(fi).folder, fns_rhd(fi).name, 0);
  dataAll = [dataAll; data];
end


bi = 3; 
flatness = flatnessAll(bi).flatness;

