function [pathout] = ZZ_winPathToLinux_v1(pathin, winDrive, linuxDrive)

%convert a windows path to linux
pathout = pathin; 
pathout(strfind(pathout,'\'))='/';
% replace the drive letter
pathout = strrep(pathout, winDrive, linuxDrive);

end

