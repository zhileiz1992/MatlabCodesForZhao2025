function [pathout] = ZZ_linuxPathToWin_v1(pathin, winDrive, linuxDrive)

%convert a windows path to linux
pathout = pathin; 
pathout(strfind(pathout,'/'))='\';
% replace the drive letter
pathout = strrep(pathout, linuxDrive, winDrive);

end

