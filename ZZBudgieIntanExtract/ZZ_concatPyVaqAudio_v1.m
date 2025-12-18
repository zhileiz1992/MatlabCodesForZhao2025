% a script to concatenate pyVAQ audio files
fd_audio =  '/mnt/z4/zz367/EphysMONAO/PyVAQData/Audio/2023-09-18';
S = dir(fullfile(fd_audio,'*.wav'));
S = S(~[S.isdir]);
% sort audio files by date
[~,idx] = sort([S.datenum]);
fns_audio = S(idx);

channelThis = 3;   % what channel to concatenate
startIdx = 119; 
endIdx = 304; 

signalTemp = [];
for ii=1:length(fns_audio)
  thisIdx = strsplit(fns_audio(ii).name, '_');
  thisIdx = strrep(thisIdx{end}, '.wav', '');
  thisIdx = str2num(thisIdx);
  if (startIdx<=thisIdx) && (thisIdx<=endIdx)
    [signal, fs] = audioread(fullfile(fns_audio(ii).folder, fns_audio(ii).name));
    signalTemp = vertcat(signalTemp, signal(:, channelThis));
  end
end

fn_save = fullfile('/mnt/z4/zz367/EphysMONAO', 'temp.wav');
audiowrite(fn_save,signalTemp, fs);
