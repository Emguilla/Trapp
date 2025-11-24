function frames2gif(frames,filename,play_dir,delay_time)
%==================================================================================================================================%
% frames2gif.m: Convert an array of MatLab frames into a GIF (v0.2)
%==================================================================================================================================%
% Version history:
%   version 0.1 (28/08/2025) - Creation
%       author: EYG
%   version 0.2 (02/09/2025) - Removal of the backward option and correction of the playing direction (it was previously going 
%       author: EYG             forward-backward-forward-forward for the infinite option)
%==================================================================================================================================%
% args:
%   frames:     Array of frames
%   filename:   Specify the filename of the GIF file 
%                   (optional argument, default: 'ReactionRendering')
%   play_dir:   Playing direction of the GIF. Possible choices are 'Fwd', 'Bwd', 'Round', 'Infinite' 
%                   (optional argument, default: 'Infinite')
%   delay_time: time delay between each frame in the GIF file
%                   (optional argument, default: 0.05)
%==================================================================================================================================%
% Check existence of optional arguments
if ~exist('filename','var')
    filename='ReactionRendering';
end
if ~exist('play_dir','var')
    play_dir='infinite';
end
if ~exist('delay_time','var')
    delay_time=0.05;
end

% Convert frames to indexed images
for p=1:length(frames)
    [img{p},c_map{p}]=rgb2ind(frame2im(frames(p)),256);
end

% set count of GIF repetition
if strcmpi(play_dir,'Infinite')
    gifcount=inf;
else
    gifcount=1;
end

% write indexed images to GIF file
imwrite(img{1},c_map{1},[filename,'.gif'],'gif','LoopCount',gifcount,'DelayTime',delay_time)
for p=[linspace(2,length(frames),length(frames)-1) linspace(length(frames)-1,2,length(frames)-2)]
    imwrite(img{p},c_map{p},[filename,'.gif'],'gif','WriteMode','append','DelayTime',delay_time)
end