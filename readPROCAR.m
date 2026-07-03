clear
close all
clc

filename='PROCAR';
path='./';

fid=fopen([path,filename],'r');
count=0;

% Header
a=fgetl(fid);
% count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
END_OF_FILE=false;
p=0;
while ~END_OF_FILE
    % nkpts, nbands and natoms informations
    a=fgetl(fid);
%    count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
    if ~ischar(a)
        END_OF_FILE=true;
    else
        p=p+1;
        Data=textscan(a,'%s %s %s %d %s %s %s %d %s %s %s %d',1,'commentStyle','%');
        nkpts=Data{4}(1);
        nbands=Data{8}(1);
        natoms=Data{12}(1);
        
        % Blank
        a=fgetl(fid);
%        count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
        
        % K-point
        for q=1:nkpts
            if q>1
                % Blank
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
            end
            a=fgetl(fid);
%            count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
            Data=textscan(a,'%s %d %s %f %f %f %s %s %f',1,'commentStyle','%');
            k{p}(q,:)=[Data{4}(1) Data{5}(1) Data{6}(1)];
            w{p}(q)=Data{9}(1);
            
            % Blank
            a=fgetl(fid);
%            count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
            
            % Bands
            for r=1:nbands
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
                Data=textscan(a,'%s %d %s %s %f %s %s %f',1,'commentStyle','%');
                E{p}(q,r)=Data{5}(1);
                Occ{p}(q,r)=Data{8}(1);
    
                % Blank
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
    
                % Band header
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
                
                % ions
                for s=1:natoms
                    a=fgetl(fid);
%                    count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
%                                   ion  s py pz px dxy dyz dz2 dxz x2-y2 tot
                    Data=textscan(a,'%d %f %f %f %f  %f  %f  %f  %f    %f  %f',1,'commentStyle','%');
                    at_orb_s{p,s}(q,r)=Data{2};
                    at_orb_py{p,s}(q,r)=Data{3};
                    at_orb_pz{p,s}(q,r)=Data{4};
                    at_orb_px{p,s}(q,r)=Data{5};
                    at_orb_dxy{p,s}(q,r)=Data{6};
                    at_orb_dyz{p,s}(q,r)=Data{7};
                    at_orb_dz2{p,s}(q,r)=Data{8};
                    at_orb_dxz{p,s}(q,r)=Data{9};
                    at_orb_x2_y2{p,s}(q,r)=Data{10};
                    at_orb_tot{p,s}(q,r)=Data{11};
                end
    
                % Total system
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
                Data=textscan(a,'%s %f %f %f %f  %f  %f  %f  %f    %f  %f',1,'commentStyle','%');
                mol_orb_s{p}(q,r)=Data{2};
                mol_orb_py{p}(q,r)=Data{3};
                mol_orb_pz{p}(q,r)=Data{4};
                mol_orb_px{p}(q,r)=Data{5};
                mol_orb_dxy{p}(q,r)=Data{6};
                mol_orb_dyz{p}(q,r)=Data{7};
                mol_orb_dz2{p}(q,r)=Data{8};
                mol_orb_dxz{p}(q,r)=Data{9};
                mol_orb_x2_y2{p}(q,r)=Data{10};
                mol_orb_tot{p}(q,r)=Data{11};
                % Blank
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
            end
        end
    end
end

fclose(fid);