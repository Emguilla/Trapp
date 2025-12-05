%function BS=readEIGENVAL(mode,varargin)
% read optional arguments
clear
close all
clc

ispath=false;
filename='EIGENVAL';
mode='bandstructure';
path='./';
% if exist('varargin')
%     for p=1:2:length(varargin)
%         switch varargin{p}
%             case 'path'
%                 ispath=true;
%                 path=varargin{p+1};
%             case 'filename'
%                 filename=varargin{p+1};
%         end
%     end
% end
% if ispath
%     if ~strcmpi(path(end),'/')&&~strcmpi(path(end),'\')
%         path=[path,'/'];
%     end
% end

EFermi=NaN;

% Face-centered cubic
%! Gamma-point %! X-point %! W-point %! L-point %! K-point %! U-point
% HS=[0 0 0; 3/8 3/8 3/4; 1/2 1/2 1/2; 5/8  1/4 5/8; 1/2 1/4 3/4; 1/2 0 1/2];
% HS_name={'\Gamma','K','L','U','W','X'};
% Hexagonal
%! Gamma-point %! L-point %! M-point
HS=[0 0 0;1/3 1/3 0;1/2 0 0];
HS_name={'\Gamma','K','M'};

if strcmpi(mode,'bandstructure')
    fid=fopen([path,filename],'r');
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    Data=textscan(fid,'%d',3,'commentStyle','%');
    nbands=Data{1}(3);
    nkpts=Data{1}(2);
    for p=1:nkpts
        fgetl(fid);
        Data=textscan(fid,'%f',4,'commentStyle','%');
        k(p,:)=Data{1}(1:3);
        for q=1:nbands
            Data=textscan(fid,'%d %f %f',1,'commentStyle','%');
            E(p,q)=Data{2};
            Occ(p,q)=Data{3};
        end
        EFermi=max([max(E(p,:).*(Occ(p,:)>0.5)) EFermi]); % def foireuse
    end
end

idx_rm=[];
HS_str={};
HS_idx=[];
is_HS=zeros(nkpts,1);
for p=1:nkpts
    for q=1:length(HS(:,1))
        if norm(k(p,:)-HS(q,:))<1e-3
            is_HS(p)=true;
            HS_str={HS_str{:},HS_name{q}};
        end
    end
end

HS_discontinuity=zeros(length(HS_str)/2+1);
HS_label{1}=HS_str{1};
for p=2:2:length(HS_str)-1
    if strcmpi(HS_str{p+1},HS_str{p})
        HS_label{p/2+1}=HS_str{p};
    else
        HS_label{p/2+1}=[HS_str{p},'/',HS_str{p+1}];
        HS_discontinuity(p/2+1)=true;
    end
end
HS_label{end+1}=HS_str{end};

BZ_int=reshape(find(is_HS),2,length(find(is_HS))/2)';
n_int=length(BZ_int(:,1));
for p=1:n_int
    for q=BZ_int(p,1):BZ_int(p,2)
        dk_int{p}(q-BZ_int(p,1)+1)=norm(k(q,:)-k(BZ_int(p,1),:));
    end
    deltak(p)=dk_int{p}(end);
    dk_int{p}=dk_int{p};
end
kx=[0 dk_int{1}(2:end)];
for p=2:n_int
    kx=[kx sum(deltak(1:p-1))+dk_int{p}(2:end)];
end
idx_rm=find(is_HS(2:end).*is_HS(1:end-1));
k(idx_rm,:)=[];
E(idx_rm,:)=[];
Occ(idx_rm,:)=[];
is_HS(idx_rm)=[];


plot(kx,E-EFermi,'k')
xlim([min(kx) max(kx)])
yl=ylim;
hold on
kx_HS=kx(find(is_HS));
for p=2:length(kx_HS)-1
    if HS_discontinuity(p)
        h=plot([kx_HS(p) kx_HS(p)],[yl(1) yl(2)],'--','color',[1 1 1]*0.75);
        uistack(h,'bottom');
    else
        h=plot([kx_HS(p) kx_HS(p)],[yl(1) yl(2)],'color',[1 1 1]*0.75);
        uistack(h,'bottom');
    end
end
xticks(kx_HS)
xticklabels(HS_label);
set(gca,'fontsize',12,'fontname','cambria math')

