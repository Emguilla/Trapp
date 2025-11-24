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
HS=[0.00000  0.00000  0.00000;%! Gamma-point
0.50000  0.00000  0.50000;%! X-point
0.50000  0.25000  0.75000;%! W-point
0.50000  0.50000  0.50000;%! L-point
0.37500  0.37500  0.75000;%! K-point
0.62500  0.25000  0.62500];%! U-point
HS_name={'\Gamma','X','W','L','K','U'};
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
        EFermi=max([max(E(p,:).*Occ(p,:)) EFermi]);
    end
end

idx_rm=[];
HS_str={};
HS_idx=[];
is_HS=zeros(nkpts,1);
for p=1:nkpts
    for q=1:length(HS(:,1))
        if norm(k(p,:)-HS(q,:))<2*eps
            is_HS(p)=true;
            HS_str=[HS_str,HS_name{q}];
            if p~=nkpts
                if norm(k(p+1,:)-k(p,:))<2*eps
                    idx_rm=[idx_rm p];
                    HS_str(end)=[];
                end
            end
        end
    end
end
k(idx_rm,:)=[];
E(idx_rm,:)=[];
Occ(idx_rm,:)=[];
is_HS(idx_rm)=[];
HS_idx=find(is_HS);
idx_break=find(is_HS(1:end-1).*is_HS(2:end));
nkx=1;
kx(1)=0;
for p=1:length(HS_idx)-1
    deltak(p)=norm(k(HS_idx(p+1),:)-k(HS_idx(p),:));
    for q=HS_idx(p)+1:HS_idx(p+1)
        nkx=nkx+1;
        kx(nkx)=kx(HS_idx(p))+(norm(k(q,:)-k(HS_idx(p),:)))*deltak(p);
    end
end
% 
% if N_ZBI_HSkpts<2
%     N_ZBI_HSkpts=Nkpoints;
%     ZBI_HSkpts_POS=1:Nkpoints;
%     ZBI_HSkpts_VAL=k;
% end
% 
% for p=1:N_ZBI_HSkpts-1
%     Nk_segments(p)=ZBI_HSkpts_POS(p+1)-ZBI_HSkpts_POS(p);
% end
% 
% intervalles=zeros(length(ZBI_HSkpts_VAL(:,1))-1,1);
% sum_intervalles=zeros(length(ZBI_HSkpts_VAL(:,1)),1);
% for p=1:length(ZBI_HSkpts_VAL(:,1))-1
%     intervalles(p)=norm(ZBI_HSkpts_VAL(p+1,:)-ZBI_HSkpts_VAL(p,:));
%     sum_intervalles(p+1)=sum(intervalles(1:p));
% end
% norm_intervalles=intervalles/sum(intervalles);
% sum_intervalles=sum_intervalles/sum(intervalles);
% 
% i=0;
% for p=1:length(ZBI_HSkpts_VAL(:,1))-1
%     for q=1:Nk_segments(p)
%         i=i+1;
%         kx=sum_intervalles(p)+(norm(ZBI_HSkpts_VAL(p,:)-k(i,:))/norm(ZBI_HSkpts_VAL(p+1,:)-ZBI_HSkpts_VAL(p,:)))*norm_intervalles(p);
%         kpts2x(i)=kx;
%     end
% end
% 
% %figure('Position', [200 120 150 425])
% for p=1:nbands
%     plot(kpts2x,E(:,p),'k')
% end
% xticklabels(ZBI_HSkpts_NAME)
% 
% hold on
% for p=2:N_ZBI_HSkpts-1
%     plot(sum_intervalles(p)*ones(2,1)-1,[min(min(E))-50 max(max(E))+50],'k-.','color',[0.75 0.75 0.75])
% end
% ylim([min(min(E)) max(max(E))])
% ylim([-25 20])
% set(gca,'fontsize',12,'fontname','cambria math')
