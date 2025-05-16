function [r_moulins] = runoff_moulins(day,runoff_2011_nevis,pp,gg)
% 23 Jan 2017 LAS
day = round(day);       % round time to nearest interger day
% mm w.e./day  to m/s
runoff_distributed = ((runoff_2011_nevis(day,:))/(1000*86400))'; % runoff on day of interest (m/sec)
runoff_repmat = repmat(runoff_distributed,1,size(pp,1))'; %(m/sec)
runoff_points = pp.*runoff_repmat; % runoff for points in catchment on day of interest (m/sec)
r_moulins = (sum(runoff_points,2)).*gg.*gg; %sum runoff for points in each moulin (m3/sec)
% catchment_size = sum(pp,2);
end
