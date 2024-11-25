% oo.yreflective = 1; oo.xreflective = 1;
oo.xperiodic = 0; oo.yperiodic = 1;
gg = nevis_grid_working2(4,4,1,0,1,0,1,oo);
gg = nevis_mask(gg,[]);
nevis_plot_grid(gg);