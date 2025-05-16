function [CLim] = splitcmapn(ax,cmaps,c)
%   Set axes in ax to have colormaps in cmaps
%   ax should be a row vector of axes handles, cmaps should be a struct with entires desired
%   colormaps for each axes, c contains handles of colorbars [optional]
%   
%   The current caxis for each axis must be large enough to encompass all
%   data within that axis [ otherwise some will be spread out onto the
%   other colormap ]

    n = length(ax);
    CLim = zeros(n,2);
    BeginSlot = zeros(n,1);
    EndSlot = zeros(n,1);
    
    % collate colormaps
    cmap = [];
    for i = 1:n
    cmap = [cmap;cmaps(i).cmap];
    end
    colormap(cmap);

    % change caxis for each axes
    for i = 1:n,
    CmLength   = length(colormap);   % Colormap length
    CLim(i,:) = get(ax(i),'CLim');
    if i==1,
    BeginSlot(i) = 1;                  % Beginning slot
    EndSlot(i)   = length(cmaps(i).cmap);      % Ending slot
    else
    BeginSlot(i) = EndSlot(i-1) + 1; 
    EndSlot(i)   = BeginSlot(i) + length(cmaps(i).cmap) - 1;
    end
    end
    for i = 1:n,
        set(ax(i),'CLim',newclim(BeginSlot(i),EndSlot(i),CLim(i,1),CLim(i,2),CmLength))
    end
    if nargin>2
    for i = 1:n
        %use number of existing Y Ticks to decide which orientation color bar is [ is there an easier way to tell the colorbar's orientation ? ]
        tmp = get(c(i),'YTick'); if isempty(tmp), set(c(i),'xlim',CLim(i,:)); else set(c(i),'ylim',CLim(i,:)); end  
    end
    end

end

function CLim = newclim(BeginSlot,EndSlot,CDmin,CDmax,CmLength)
   %     Convert slot number and range
   %     to percent of colormap
   PBeginSlot    = (BeginSlot - 1) / (CmLength - 1);
   PEndSlot      = (EndSlot - 1) / (CmLength - 1);
   PCmRange      = PEndSlot - PBeginSlot;
   %     Determine range and min and max 
   %     of new CLim values
   DataRange     = CDmax - CDmin;
   ClimRange     = DataRange / PCmRange;
   NewCmin       = CDmin - (PBeginSlot * ClimRange);
   NewCmax       = CDmax + (1 - PEndSlot) * ClimRange;
   CLim          = [NewCmin,NewCmax];
end
    