function [rwdarms] = rwds(filen)
% Function to return the rewarded arms for each animal;
% place-to-go project
% adamhsugi@gmail.com

%L M H

rat = filen(2:3);

if strcmp(rat,'09')
        rwdarms = [3 ,6 ,8];
elseif strcmp(rat,'08')
        rwdarms = [2 ,5 ,7];
elseif strcmp(rat,'11')
        rwdarms = [3 ,8 ,5];
elseif strcmp(rat,'19')
        rwdarms = [3 ,8 ,5];
elseif strcmp(rat,'04')
        rwdarms = [2 ,5 ,7];
end 