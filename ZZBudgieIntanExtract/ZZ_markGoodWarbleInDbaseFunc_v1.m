% A function to select subsets of files in a dbase
% focus on the warble episodes of focal bird
% modified to work with Brian's new gui
function [dbase_new, dbase_old] = ZZ_markGoodWarbleInDbaseFunc_v1(fn_dbase, fn_select)
%   fn_dbase = dir(fullfile(fd_base_dbase, expID, data_dates{di}, birdID, 'warble', sprintf('%s.*.warble.segFlatness.dbase.mat', birdID)));
%   load(fullfile(fn_dbase.folder, fn_dbase.name));
  load(fn_dbase);
  % what warble episodes to select
%   fn_select = dir(fullfile(fd_base_dbase, expID, data_dates{di}, birdID, 'warble', sprintf('%s.*.warble.select.txt', birdID))); 
  lines = readlines(fn_select);
  to_select = [];
  for li=1:length(lines)
    ln_rep = strrep(lines{li}, ',', ' ');
    array = eval(['[', ln_rep, ']']);
    to_select = [to_select array];
  end
  % determine those warble episode corresponds to what dbase files
  fns_all = {dbase.SoundFiles(:).name};
  select_idx = [];
  for fi=1:length(fns_all)
    % get the warble id
    wid = strsplit(fns_all{fi}, '_');
    wid = str2num(wid{2});
    if ismember(wid, to_select)
      select_idx = [select_idx fi];
    end
  end
  % mark the bGood property in dbase
  dbase_old = dbase;
  dbase_new = dbase;
  % which property id bGood property
  good_idx = find(strcmp(dbase.PropertyNames, 'bGood'));
  a = dbase.Properties;
  a(select_idx, good_idx) = true;
  dbase_new.Properties = a;

end

  
    
  
  