function [ni_g] = nevis_gps_array(xg,yg,gg,oo)
    ni_g = zeros(1,length(xg)*length(yg)); % Initialize GPS array indices
    for i = 1:length(xg)
        for j = 1:length(yg)
            % Ensure xg and yg are column vectors
            vxg = xg(i); % 
            vyg = yg(j); %

            % find the closest index to each GPS point
            [~,tmp] = min((reshape(gg.nx,[],1)-vxg).^2 + ...
                        (reshape(gg.ny,[],1)-vyg).^2);
            ni_g(i+(j-1)*length(xg)) = tmp; % store the index
        end
    end
    ni_g = unique(ni_g, 'stable'); % ensure unique indices
end