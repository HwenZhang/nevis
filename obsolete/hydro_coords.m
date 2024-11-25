function [nx,ny,ex,ey,fx,fy,cx,cy] = hydro_coords(nI,nJ,xl,xr,yb,yt,opts)
% coordinates of grid points
% [ grid is xl<x<xr, yb<y<yt ]
%
% 9 May 2013 : added option for reflective boundary conditions in y

%one of the following methods for choosing boundaries should be chosen :
if exist('opts','var') && isfield(opts,'method') method = opts.method; else
%   method = 'lrtb_vel';               % appropriate for prescribing velocity on all boundaries 
%   method = 'lr_stress_tb_vel';       % appropriate for prescribing velocity on top and bottom boundaries, stress on left and right 
%   method = 'r_stress_ltb_vel';       % appropriate for prescribing velocity on top, bottom and left boundaries, stress on right
   method = 'lrtb_per';           % appropriate for periodic top, bottom, left and right boundaries [ may need to fiddle periodic edges ]
%   method = 'lr_per_tb_vel';          % appropriate for prescribing velocity on top and bottom boundaries, periodic left and right 
%   method = 'r_stress_l_vel_tb_per';  % appropriate for periodic top and bottom boundaries, velocity on left, stress on right 
%   method = 'lr_stress_tb_per';       % appropriate for periodic top and bottom boundaries, stress on left and right [ needs checking ]
%   method = 'mask_stress';            % use some condition to mask the domain and prescribe stress on all boundaries
%   method = 'mask_vel';               % use some condition to mask the domain and prescribe velocity on all boundaries
%   method = 'r_stress_l_vel_tb_ref';  % appropriate for reflective top and bottom boundaries,
%   velocity on left, stress on right 
end

%     % [ this will make edges on boundaries ]
%     ex = (0:nI)/nI;
%     fy = (0:nJ)/nJ;
%     nx = (ex(1:end-1)+ex(2:end))/2;
%     ny = (fy(1:end-1)+fy(2:end))/2;
%     ey = ny;
%     fx = nx;
%     cx = ex;
%     cy = fy;

if strcmp(method,'lrtb_vel') || strcmp(method,'lr_stress_tb_vel') || strcmp(method,'r_stress_ltb_vel') || strcmp(method,'mask_stress') || strcmp(method,'mask_vel') || strcmp(method,'r_stress_l_vel_tb_ref')
    % [ this will make nodes on boundaries ]
    ex = (-1/2+(0:nI))/(nI-1);
    fy = (-1/2+(0:nJ))/(nJ-1);
    nx = (ex(1:end-1)+ex(2:end))/2;
    ny = (fy(1:end-1)+fy(2:end))/2;
    ey = ny;
    fx = nx;
    cx = ex;
    cy = fy;

elseif strcmp(method,'r_stress_l_vel_tb_per') || strcmp(method,'lr_stress_tb_per')
    % [ this will make nodes on boundaries at left and right and bottom, for periodic bcs in y ]
    ex = (-1/2+(0:nI))/(nI-1);
    fy = (-1/2+(0:nJ))/(nJ);
    nx = (ex(1:end-1)+ex(2:end))/2;
    ny = (fy(1:end-1)+fy(2:end))/2;
    ey = ny;
    fx = nx;
    cx = ex;
    cy = fy;
    
elseif strcmp(method,'lrtb_per') || strcmp(method,'lr_per_tb_vel')
    % [ this will make nodes on boundaries at left and bottom, for periodic bcs in x and y ]
    ex = (-1/2+(0:nI))/(nI);
    fy = (-1/2+(0:nJ))/(nJ);
    nx = (ex(1:end-1)+ex(2:end))/2;
    ny = (fy(1:end-1)+fy(2:end))/2;
    ey = ny;
    fx = nx;
    cx = ex;
    cy = fy;
    
else disp('hydro_coords : Unsupported boundary conditions');
    
end    
    
    % scale to size
    nx = xl+(xr-xl)*nx;
    ny = yb+(yt-yb)*ny;
    ex = xl+(xr-xl)*ex;
    ey = yb+(yt-yb)*ey;
    fx = xl+(xr-xl)*fx;
    fy = yb+(yt-yb)*fy;
    cx = xl+(xr-xl)*cx;
    cy = yb+(yt-yb)*cy;

    [ny,nx] = meshgrid(ny,nx);
    [ey,ex] = meshgrid(ey,ex);
    [fy,fx] = meshgrid(fy,fx);
    [cy,cx] = meshgrid(cy,cx);
    
end