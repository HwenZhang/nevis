function transf_fn = get_transform(proj1,proj2)
%  transf_fn = get_transform(proj1, proj2);
%
% Returns transformation function to transform coordinates in proj1 to coordinates in
% proj2.  The returned function has the following signature:
% [x2,y2] = transf_fn(x1,y1);
% which can be called with scalars, vectors or matrices.
%
% Can do:
% - polar stereo to polar stereo
% - latlong to polar stereo and vice versa
% 
% Input proj requires
%   proj.type   'ps' or 'latlog' or 'rot_latlong', or specific option as below
%   'ps'
%   proj.earthradius
%   proj.eccentricity 
%   proj.std_parallel 
%   proj.central_meridian
%   'rot_latlong' lat-long rotated pole http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/ch05s06.html\
%   proj.north_pole_long-180, 
%   proj.north_pole_lat
%
% 23 Sept 2014: modified from Mauro Werder's get_transform.m

% default constants:
WGS84_earthradius = 6378137.0;
WGS84_eccentricity = 0.08181919;

% specific projection types 
for proj = [proj1 proj2]
if strcmp(proj.type,'bamber') || strcmp(proj.type,'searise')
% polar stereo: bamber-like
proj.type = 'ps';
proj.units = 'm';
proj.earthradius = WGS84_earthradius;
proj.eccentricity =  WGS84_eccentricity;
proj.std_parallel = 71;
proj.origin = 90;
proj.central_meridian = -39;
elseif strcmp(proj.type,'gimp') || strcmp(proj.type,'joughin')
% polar stereo: gimp-like (aka ESPG 3413)
proj.type = 'ps';
proj.units = 'm';
proj.earthradius = WGS84_earthradius;
proj.eccentricity =  WGS84_eccentricity;
proj.std_parallel = 70;
proj.origin = 90;
proj.central_meridian = -45;
elseif strcmp(proj.type,'antarctica') 
% polar stereo: Antarctica-like]
proj.type = 'ps';
proj.units = 'm';
proj.earthradius = WGS84_earthradius;
proj.eccentricity =  WGS84_eccentricity;
proj.std_parallel = -70;
proj.origin = -90;
proj.central_meridian = 0;
end
end

switch proj1.type
  case 'none'
    if proj2.type=='none'
        transf_fn = @id;
        return;
    else
        error('cannot do transformation');
    end
  case 'ps' % polar stereo
    switch proj2.type
      case 'ps'
        transf_fn = ps2ps(proj1, proj2);
      case 'latlong'
        transf_fn = ps2ll(proj1, proj2);
      case 'rot_latlong'
        transf_fn = ps2rot(proj1, proj2);
      otherwise
        error('cannot do transformation');
    end
  case 'latlong' % lat long
    switch proj2.type
      case 'ps'
        transf_fn = ll2ps(proj1, proj2);
      case 'latlong'
        transf_fn = @id;
      case 'rot_latlong'
        transf_fn = ll2rot(proj1, proj2);
      otherwise
        error('cannot do transformation');
    end
  case 'rot_latlong'
    switch proj2.type
      case 'latlong'
        transf_fn = rot2ll(proj1, proj2);
      case 'ps'
        transf_fn = rot2ps(proj1, proj2);
      case 'rot_latlong'
        transf_fn = rot2rot(proj1, proj2);
      otherwise
        error('cannot do transformation');
    end
  otherwise 
    error('cannot do transformation');
end
end

%%%%%%%%%%%%%%%%%%%%
% Function factories
%
% These are needed because anonymous functions cannot return two values, and we want to
% hide the input p1, p2.

function fn = ps2ps(p1,p2)
    fn = @(x, y) ps2ps_fn(x, y, p1, p2);
end 

function fn = ll2ps(p1,p2)
    fn = @(x, y) ll2ps_fn(x, y, p1, p2);
end

function fn = ps2ll(p1,p2)
    fn = @(x, y) ps2ll_fn(x, y, p1, p2);
end

function fn = rot2ll(p1,p2)
    fn = @(x,y) rot2ll_fn(x, y, p1, p2);
end

function fn = ll2rot(p1,p2)
    fn = @(x,y) ll2rot_fn(x, y, p1, p2);
end

function fn = ps2rot(p1,p2)
    fn = @(x,y) ps2rot_fn(x, y, p1, p2);
end

function fn = rot2ps(p1,p2)
    fn = @(x,y) rot2ps_fn(x, y, p1, p2);
end

function fn = rot2rot(p1,p2)
    fn = @(x,y) rot2rot_fn(x, y, p1, p2);
end


%%%%%%%%%%%%%%%%%%%%
% Elementary transformation functions
function [x, y] = id(x, y, ~, ~)
    x=x; y=y;
end

function [x, y] = ll2ps_fn(x, y, ~, p2)
    [x, y]= polarstereo_fwd(x, y, p2.earthradius, p2.eccentricity, ...
                            p2.std_parallel, p2.central_meridian);
end

function [x, y] = ps2ll_fn(x, y, p1, ~)
    [x, y]= polarstereo_inv(x, y, p1.earthradius, p1.eccentricity, ...
                            p1.std_parallel, p1.central_meridian);
end

function [x,y] = rot2ll_fn(x, y, p1, ~)
    option = 2;
    % south pool lat/long
    SP_coord = [p1.north_pole_long-180, -p1.north_pole_lat];
    [x,y] = rotated_grid_transform(x, y, option, SP_coord);
end

function [x,y] = ll2rot_fn(x, y, ~, p2)
    option = 1;
    % south pool lat/long
    SP_coord = [p2.north_pole_long-180, -p2.north_pole_lat];
    [x,y] = rotated_grid_transform(x,y, option, SP_coord);
end

%%%%%%%%%%%%%%%%%%%%%%%55
% Composite transformation functions

function [x, y] = ps2ps_fn(x, y, p1, p2)
    [x,y] = ps2ll_fn(x, y, p1, p2);
    [x, y] = ll2ps_fn(x, y, p1, p2);
end 

function [x, y] = rot2rot_fn(x, y, p1, p2)
    [x,y] = rot2ll_fn(x, y, p1, p2);
    [x,y] = ll2rot_fn(x, y, p1, p2);
end 

function [x, y] = rot2ps_fn(x, y, p1, p2)
    [x,y] = rot2ll_fn(x, y, p1, p2);
    [x, y] = ll2ps_fn(x, y, p1, p2);
end 

function [x, y] = ps2rot_fn(x, y, p1, p2)
    [x,y] = ps2ll_fn(x, y, p1, p2);
    [x,y] = ll2rot_fn(x, y, p1, p2);
end 

