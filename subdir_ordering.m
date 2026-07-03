function [ldir,lname,lsname]=subdir_ordering(dirpath)
%==================================================================================================================================%
% subdir_ordering.m:   Recursive navigation of NEB and subNEB directories (v0.1.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (01/07/2026) - Creation
%       author: EYG
%   version 0.1.1 (02/07/2026) - Add a recursive definition of the image depth
%       contrib: EYG
%   version 0.1.2 (03/07/2026) - The first and last image in the subNEB directories are removed only if a folder named "00" is
%       contrib: EYG                present.
%==================================================================================================================================%
% args:
%   dirpath:        Location of the main directory
%==================================================================================================================================%

% Save current directory location to get back to the same point after execution, then jumps to the specified directory
curr_dir=pwd;
cd(dirpath)
image_depth=[];
% List everything in the directory
ldir=dir;
% Only keep sub-directories
ldir=ldir([ldir(:).isdir]);
% Remove everything that is not a numbered directory (i.e. we only keep the images)
clear_idx=[];
for p=1:length(ldir)
    if isempty(str2num(ldir(p).name))
        clear_idx=[clear_idx p];
    end
end
ldir(clear_idx)=[];

% Generation of the list of directories as string of numbers
for p=1:length(ldir)
    lname{p}=ldir(p).name;
    image_depth{p}=0; % Here we define the depth of the image. It is set to zero, then will be incremented with each recursive call
end
[ldir.image_depth]=image_depth{:};
% Check how many padding zeros are present
format_length=length(ldir(1).name);
format=['%0',num2str(format_length),'i'];

% Identify subNEB directories
target_dir=dir('subNEB_*');
rmdir_count=0;
if ~isempty(target_dir)
    for p=1:length(target_dir)
        clear lsubname
        % Recursive call to the function to analyse subNEB directories within subNEB directories
        [lsubdir,lsubname_tmp]=subdir_ordering(target_dir(p).name);
        % Add one to the depth of the subNEB images
        for q=1:length(lsubdir)
            lsubdir(q).image_depth=lsubdir(q).image_depth+1;
        end
        % Removal of the endpoints of the subNEB as they are the same as the one from the parent directory, only if the "00"
        % directory exists
        if str2double(lsubdir(1).name)==0
            lsubdir([1 end])=[];
            lsubname_tmp([1 end])=[];
        end
        % The most common case for a subNEB is to run it between consecutive images. However, if you ran the subNEB by covering two
        % or three intervals (like BC between 01 and 03, skipping 02; or DEF between 03 and 06, skipping 04 and 05), then the 
        % intermediate(s) image(s) is(are) removed from the list.
        if length(target_dir(p).name)==8
            idx_up(p,:)=[double(target_dir(p).name(end))-64 double(target_dir(p).name(end))-63];
            txt_prefix=target_dir(p).name(end-7:end);
        elseif length(target_dir(p).name)==9
            idx_up(p,:)=[double(target_dir(p).name(end-1))-64 double(target_dir(p).name(end))-63];
            txt_prefix=target_dir(p).name(end-8:end);
            ldir(idx_up(p,1)+1-rmdir_count)=[];
            lname(idx_up(p,1)+1-rmdir_count)=[];
            rmdir_count=rmdir_count+1;
        elseif length(target_dir(p).name)==10
            idx_up(p,:)=[double(target_dir(p).name(end-2))-64 double(target_dir(p).name(end))-63];
            txt_prefix=target_dir(p).name(end-9:end);
            ldir(idx_up(p,1)+1-rmdir_count)=[];
            lname(idx_up(p,1)+1-rmdir_count)=[];
            rmdir_count=rmdir_count+1;
            ldir(idx_up(p,1)+2-rmdir_count)=[];
            lname(idx_up(p,1)+2-rmdir_count)=[];
            rmdir_count=rmdir_count+1;
        end
        % The subNEB_X is replaced by the number of the image used as reactant of the subNEB
        prefix=[num2str(idx_up(p,1)-1,format),'/'];
        for q=1:length(lsubname_tmp)
            lsubname{q}=[prefix,lsubname_tmp{q}];
            lsubdir(q).name=[txt_prefix,'/',lsubdir(q).name];
        end
        ldir(length(ldir)+1:length(ldir)+length(lsubdir))=lsubdir;
        lname=[lname lsubname];
    end
end

% Thanks to the fact we converted the subNEB to the number of the reactants, we can simply "sort" the list of directories
[lname,idx]=sort(lname);
ldir=ldir(idx);

% For shorthand notation, the "subNEB" is removed from the name of the paths.
for p=1:length(ldir)
    lsname{p}=strrep(ldir(p).name,'subNEB_','');
end

% Relocate to initial directory
cd(curr_dir)