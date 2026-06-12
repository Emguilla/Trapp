clear
close all
clc

filename='PROCAR';
path='./';

fid=fopen([path,filename],'r');
count=0;

% Header
a=fgetl(fid);
count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
END_OF_FILE=false;
p=1;
while ~END_OF_FILE
    % nkpts, nbands and natoms informations
    a=fgetl(fid);
    count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
    if a==-1
        END_OF_FILE=true;
    else
        p=p+1;
        disp(p)
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
                E(p,q)=Data{5}(1);
                Occ(p,q)=Data{8}(1);
    
                % Blank
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
    
                % Band header
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
                
                % ions
                for p=1:natoms
                    a=fgetl(fid);
%                    count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
                end
    
                % Total system
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a])
        
                % Blank
                a=fgetl(fid);
%                count=count+1;disp([' reading #',num2str(count,'%03i'),': ',a,' BLANK'])
            end
        end
    end
end

fclose(fid);