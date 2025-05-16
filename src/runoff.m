function [r] = runoff(day,runoff_2011_nevis)
% 10 May 2016 LAS
day = round(day);       % round time to nearest interger day
% mm w.e./day  to m/s
r = ((runoff_2011_nevis(day,:))/(1000*86400))'; % runoff on day of interest (m/sec)
end