function outfile = nevis_gps_stress(casename, varargin)
% Thin wrapper around nevis_gps_stresses for a singular function entrypoint.

outfile = nevis_gps_stresses(casename, varargin{:});

end