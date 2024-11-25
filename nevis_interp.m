function out = nevis_interp(in,x_in,y_in,proj_in,x,y,proj,oo)
% interpolate quantity in defined on mesh x_in,y_in in projection proj_in, to mesh x,y in projection
% proj
%
% Example interpolate speed v defined on grid [x1,y1] in proj1 onto grid [x,y] in proj2
%  proj1.type = 'joughin'; proj2.type = 'bamber';
%  v = nevis_interp(v,x_in,y_in,proj1,x,y,proj2);
%
% 23 Sept 2014: modified from Mauro Werder's interp_geo.m [ might want to include option for
%   temporal interpolation here too ]
% 3 Nov 2015: changed default from unstructured input to structured, as
%   griddata seems slow

if nargin<8, oo = struct; end
if ~isfield(oo,'structured_input'), oo.structured_input = 1; end 
if size(x_in)~=size(in), [x_in,y_in] = meshgrid(x_in,y_in); end % meshgrid input if x_in and y_in included as vectors?

% transformation to input grid
if ~strcmp(proj_in.type,proj.type),
[x2,y2] = nevis_transform(x,y,proj,proj_in);
else
x2 = x; y2 = y;
end

% crop input grid to relevant pts
xmin = min(min(x2)); xmax = max(max(x2));
ymin = min(min(y2)); ymax = max(max(y2));

if oo.structured_input
% structured input grid
xpts = x_in(:,1)>=xmin & x_in(:,1)<=xmax;
ypts = y_in(1,:)>=ymin & y_in(1,:)<=ymax;
out = interp2(x_in(xpts,ypts)',y_in(xpts,ypts)',in(xpts,ypts)',x2',y2','linear',NaN)';
% out = interp2(x_in,y_in,in,x2,y2,'linear',NaN);
else
% unstructured input grid
pts = x_in>=xmin & x_in<=xmax & y_in>=ymin & y_in<=ymax;
out = griddata(x_in(pts),y_in(pts),in(pts),x2,y2,'linear');
end

end