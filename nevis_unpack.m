function nevis_unpack(varargin)
% e.g. nevis_unpack(gg,vv)
% uupacks the contents of structures to the current workspace
    for is = 1:length(varargin)
        s = varargin{is};
        names = fieldnames(s);
        for in = 1:length(names),
            assignin('caller',names{in},getfield(s,names{in}));
        end
    end
end