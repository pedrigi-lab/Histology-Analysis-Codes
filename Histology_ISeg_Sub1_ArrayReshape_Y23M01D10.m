% This subroutine identifies the files needed within a folder and organizes them.
%
% R. Pedrigi (Dec, 2018)

function filedirimage = Histology_ISeg_Sub1_ArrayReshape_Y23M01D10(filedir1,datadir1,mouseid,stain,filestr)

%%% For running standalone
% filestr = 'M*';
% filestr = 'Phase*';
% filestr = 'Binary*';


if strcmp(filestr,'Phase*') || strcmp(filestr,'Binary*')
    filebasename = strcat(filestr(1:end-1),'_',mouseid);
else 
    filebasename = mouseid;
end
 
% Find all raw images for the mouse stain.
if strcmp(filestr,'M*') == 1
    dirfinal = filedir1;
else
    dirfinal = datadir1;
end

listing = dir(strcat(dirfinal,filestr,'.tif'));
filediralltmp = cell(1,1);
num = 0;
for i = 1:size(listing,1)
    
    tmp = listing(i).name;             % List all files in data directory (if there are none, it is accounted for below).
    if strcmp(tmp(size(tmp,2)-2:size(tmp,2)),'mat') == 0
        num = num+1;
        filediralltmp{num,:} = strcat(dirfinal,listing(i).name);
    end
    
end


% Create strings for the different section series, in order
strtmp1 = strcat(dirfinal,filebasename,'L_',stain,'_SI_');
strtmp2 = strcat(dirfinal,filebasename,'L_',stain,'_SII_');
strtmp3 = strcat(dirfinal,filebasename,'L_',stain,'_SIII_');
strtmp4 = strcat(dirfinal,filebasename,'L_',stain,'_SIV_');
strtmp5 = strcat(dirfinal,filebasename,'L_',stain,'_SV_');
strtmp6 = strcat(dirfinal,filebasename,'L_',stain,'_SVI_');

% Reorganize the images so they are in sequence
strorder = zeros(size(filediralltmp,1),1);
filedirimage = cell(size(filediralltmp,1),1);
i = 0; k = 0;
tmp = 0;
while tmp < size(filediralltmp,1)
    
    i = i+1;
    
    strtmp0 = eval(strcat('strtmp',num2str(i)));
    
    strorder(:,i) = contains(filediralltmp,strtmp0);
    
    tmp = tmp+sum(strorder(:,i));
    
    j = k+1;
    k = k+sum(strorder(:,i));
    for m = 1:sum(strorder(:,i))
        loc = find(strorder(:,i)==1,1);
        filedirimage{(j-1)+m,1} = filediralltmp{loc+(m-1),1};
    end
    
end







