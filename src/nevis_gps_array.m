function [ni_g] = nevis_gps_array(xg,yg,gg,oo)
    ni = min(length(xg), length(yg)); % total number of grid points
    ni_g = zeros(1,ni); % Initialize GPS array indices
    for i = 1:ni
        % Ensure xg and yg are column vectors
        vxg = xg(i); % 
        vyg = yg(i); %

        % find the closest index to each GPS point
        [~,tmp] = min((reshape(gg.nx,[],1)-vxg).^2 + ...
                    (reshape(gg.ny,[],1)-vyg).^2);
        ni_g(i) = tmp; % store the index
    end
    % ni_g = unique(ni_g, 'stable'); % ensure unique indices
end