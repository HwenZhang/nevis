function is = nevis_orderboundary(x,y)
% find indices to order boundary nodes x,y in approximate tracing order
% does not always work!
% 
% 26 August 2013

    % options
    large_d_fac = 5; % 5
    d_fac = 1.5; % 1.5

    x = reshape(x,[],1); 
    y = reshape(y,[],1);
    xx = x*x'.^0-x.^0*x';
    yy = y*y'.^0-y.^0*y';
    dd = (xx.^2+yy.^2).^(1/2);
    
    is = NaN*ones(1,length(x)); 
    ds = NaN*ones(1,length(x));
    ee = []; % dead ends
    i1 = 1; is(1) = 1; 
    rr = setdiff(1:length(x),is(1));
    while length(rr)>0
        i2 = is(i1);
        [d,tmp] = min(dd(i2,rr)); ds(i1) = d; i3 = rr(tmp);
        if d>large_d_fac*mean(ds(1:i1)), is(i1+1:end) = is(i1); return; end % stop if too large a step taken since that suggests pts that were earlier left out
        if d>d_fac*mean(ds(1:i1)), is = [is is(end)]; ds = [ds ds(end)]; 
            ee = union(ee,i2);
        rr2 = setdiff(1:length(x),ee); [d,tmp] = min(dd(i2,rr2)); ds(i1) = d; i3 = rr2(tmp); end % allow retracing steps if all pts are too far away

        is(i1+1) = i3;
        i1 = i1+1;
        rr = setdiff(1:length(x),is(1:i1));
%         plot(x,y,'.',x(is(1:i1)),y(is(1:i1)),'o'); shg;
    end
end