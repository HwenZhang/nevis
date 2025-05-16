function i_out = int2four(i)    
    i1 = floor(i/1000); i2 = floor((i-1000*i1)/100); i3 = floor((i-1000*i1-100*i2)/10); i4 = i-1000*i1-100*i2-10*i3;
    i_out = [num2str(i1),num2str(i2),num2str(i3),num2str(i4)];
end
