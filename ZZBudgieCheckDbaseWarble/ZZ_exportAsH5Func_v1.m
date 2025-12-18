function fd_h5 = ZZ_exportAsH5Func_v1(sound_struct, spec_field, fd_h5, spec_size, prefix)
%ZZ_EXPORTASH5FUNC_V1 Summary of this function goes here
% export the spectrograms in sound_struct as H5 files for network traing
% loop through all syllables, export as h5 file
if exist(fd_h5)
  rmdir(fd_h5, 's')
end
mkdir(fd_h5);

parfor si=1:length(sound_struct)
  fn = fullfile(fd_h5, sprintf('%s_%07d.h5', prefix, si));
  h5create(fn, '/spec', spec_size);
  h5writeatt(fn,'/spec','budgieID', prefix);
  % interpolate the spectrogram to desired size
  A = sound_struct(si).(spec_field);
  [n, m] = size(A); 
  [x, y] = meshgrid(1:m, 1:n);
  [xi, yi] = meshgrid(linspace(1, m, spec_size(2)), linspace(1, n, spec_size(1)));
  Ai = interp2(x, y, A, xi, yi, 'linear');
  h5write(fn,'/spec',Ai);
end
end

