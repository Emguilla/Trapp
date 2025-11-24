function x=coordinate_mapping(XDATCARs,parametric_coordinates,projected_coordinates,projected_resolution)
%==================================================================================================================================%
% coordinate_mapping.m: computation of the reaction coordinates of the POSCAR structure matrix corresponding to an ensemble of 
%                       XDATCAR files (v0.1)
%==================================================================================================================================%
% Version history:
%   version 0.1 (11/09/2025) - Creation
%       author: EYG
%==================================================================================================================================%
% args:
%   XDATCARs:   Matrix of POSCAR structures 
%                   (size along dimension 1: number of images)
%                   (size along dimension 2: number of iterations)
%==================================================================================================================================%
n_images=length(XDATCARs(:,1))-2;
n_iter=length(XDATCARs(1,:));
% Parametric mapping
if parametric_coordinates
    % for parametric coordinates, only the distance between successive points is computed, and the coordinates of the i-th images 
    % is the sum of all the successive distance up to that image
    ds=zeros(n_images+2,n_iter);
    for q=1:n_iter
        for p=2:n_images+2
            if p==2
                ds(p,q)=ds_POSCAR(XDATCARs(p,q),XDATCARs(p-1,end));
            elseif p==n_images+2
                ds(p,q)=ds_POSCAR(XDATCARs(p,end),XDATCARs(p-1,q));
            else
                ds(p,q)=ds_POSCAR(XDATCARs(p,q),XDATCARs(p-1,q));
            end
        end
        for p=1:n_images+2
            x(p,q)=sum(ds(1:p,q))./sum(ds(:,q));
        end
    end
% Projection mapping
elseif projected_coordinates
    m=0:projected_resolution:1;
    n_interp=length(m);
    % It is highly preferable to put a POSCAR_mid file in the folder or its parent to explicitily choose the intermediate point. If
    % not, a rough analysis of the images is performed to assess which one is farthest away from both the reactants and products.
    if ~exist('POSCAR_mid','file')&&~exist('../POSCAR_mid','file')
        for p=1:n_images+2
            if p~=1&&p~=n_images+2
                Fwd_ds(p)=ds_POSCAR(XDATCARs(1,1),XDATCARs(p,end));
                Bwd_ds(p)=ds_POSCAR(XDATCARs(end,1),XDATCARs(p,end));
            else
                Fwd_ds(p)=ds_POSCAR(XDATCARs(1,1),XDATCARs(p,1));
                Bwd_ds(p)=ds_POSCAR(XDATCARs(end,1),XDATCARs(p,1));
            end
        end
        [~,idx_max_ds]=max(Bwd_ds+Fwd_ds);
        pseudoTS=XDATCARs(idx_max_ds,end);
        warning(sprintf(['It is proferable to provide a guess of the intermediate configuration\n' ...
            'geometry in order to find the projected coordinate\n' ...
            'Proceed with caution as the default will probably not be adequate !']))
    else
        if ~exist('POSCAR_mid','file')
            pseudoTS=readPOSCAR('../POSCAR_mid');
        else
            pseudoTS=readPOSCAR('POSCAR_mid');
        end
        % Find to what index correspond the closest image to the POSCAR_mid configuration
        for p=1:n_images+2
            ds_pTS(p)=ds_POSCAR(pseudoTS,XDATCARs(p,end));
        end
        [~,idx_max_ds]=min(ds_pTS);
    end
    % get distance and vector to get from reactant to TS, then from TS to products
    [ds_A,ds_A_vec,ds_A_xred]=ds_POSCAR(XDATCARs(1,1),pseudoTS);
    [ds_B,ds_B_vec,ds_B_xred]=ds_POSCAR(pseudoTS,XDATCARs(end,1));
    ds_tot=ds_A+ds_B;
    % Interpolation along the line between reactants and TS, then between TS and products
    for r=1:n_interp
        m_interp_A(r)=XDATCARs(1,1);
        m_interp_A(r).positions=m_interp_A(r).positions+m(r)*ds_A_vec;
        m_interp_A(r).xred=m_interp_A(r).xred+m(r)*ds_A_xred;
        m_interp_B(r)=pseudoTS;
        m_interp_B(r).positions=m_interp_B(r).positions+m(r)*ds_B_vec;
        m_interp_B(r).xred=m_interp_B(r).xred+m(r)*ds_B_xred;
    end
    x=zeros(n_images+2,n_iter);
    x(end,:)=1;
    for q=1:n_iter
        % for each iteration, calculation of the distance between the "soon-to-be-TS" and reactants and products
        ds_g=ds_POSCAR(XDATCARs(1,1),XDATCARs(idx_max_ds,q));
        ds_r=ds_POSCAR(XDATCARs(end,1),XDATCARs(idx_max_ds,q));
        for p=2:n_images+1
            if p<idx_max_ds % i.e. left-side of the reaction
                % Get distances between the current state of an image and the interpolated points between reactant and TS
                for r=1:n_interp
                    ds(r)=ds_POSCAR(XDATCARs(p,q),m_interp_A(r));
                end
                % The index of the minimum of the distances between the current state of the image and the interpolated
                % configurations indicates the relative progression (between 0 and n_interp-1) along the reaction coordinates
                [~,idx_min_ds]=min(ds);
                x(p,q)=ds_g*(idx_min_ds-1)/((n_interp-1)*(ds_g+ds_r));
            elseif p>idx_max_ds % i.e. right side of the reaction
                % Ditto case p<idx_max_ds
                for r=1:n_interp
                    ds(r)=ds_POSCAR(XDATCARs(p,q),m_interp_B(r));
                end
                [~,idx_min_ds]=min(ds);
                x(p,q)=ds_g/(ds_g+ds_r)+ds_r*(idx_min_ds-1)/((n_interp-1)*(ds_g+ds_r));
            else
                % for the "soon-to-be-TS" images, the reaction coordinates is the ratio between the distance to reactants over the
                % total distance reactants/TS/products
                x(p,q)=ds_g/(ds_g+ds_r);
            end
        end
    end
end

