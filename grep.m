function output_line_txt=grep(filename,string,varargin)
%==================================================================================================================================%
% grep.m: Extraction of a specific pattern in a file (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (20/08/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   filename:   path to file including its name
%   string:     pattern to search for in the file
%   opt. args:  'fwd', followed by an integer to request a certain amount of line after to the appearance of the pattern
%               'bwd', followed by an integer to request a certain amount of line prior to the appearance of the pattern
%               'prt', followed by a character string is an option to write the output to the file 
%               'vis', followed by true or false (default) is a verbose option to display the output on screen
%==================================================================================================================================%
% Initialisation of the default parameters
output_line_txt={};
fwd=0;
tbprt=false;
vis=false;
% read optional arguments
if exist('varargin')
    for p=1:2:length(varargin)
        switch varargin{p}
            case 'fwd'
                fwd=varargin{p+1};
            case 'bwd'
                bwd=varargin{p+1};
                disp('not yet implemented')
            case 'prt'
                prt=varargin{p+1};
                tbprt=true;
            case 'vis'
                vis=varargin{p+1};
        end
    end
end

% opening of the file in which the search must be done
fid=fopen(filename);
k=1;
line_txt=fgetl(fid);
% for each line of the file, search for a match with the pattern. If successful, the line containing the match is stored in a cell
% array
while ischar(line_txt)
    if contains(line_txt,string)
        output_line_txt{k}=line_txt;
        if vis
            disp(line_txt)
        end
        k=k+1;
        if fwd~=0
            p=1;
            while p<=fwd
                line_txt=fgetl(fid);
                output_line_txt{k}=line_txt;
                if vis
                    disp(line_txt)
                end
                k=k+1;
                if contains(line_txt,string)
                    p=p-1;
                end
                p=p+1;
            end
        end
    end
    line_txt=fgetl(fid);
end
fclose(fid);
% If relevant, the output is written to a file
if tbprt
    fid=fopen(prt,'w');
    fprintf(fid,'%s \n',output_line_txt{:});
    fclose(fid);
end
