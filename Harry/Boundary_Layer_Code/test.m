% animateThinFilm.m
% Call Composite_Solution to create a thin film evolution animation

% Physical and numerical parameters
h0 = 1;             % thickness scale (m)
L0 = 1000;          % length scale (m)
t0 = 3600;          % time scale (s)
rhow = 1000;        % density of water (kg/m^3)
g = 9.81;           % gravitational acceleration (m/s^2)
mu = 1e1;           % dynamic viscosity of water (Pa.s)
theta = sin(0.1);   % slope angle (radians)


% —— Parameter settings —— 
Ll    = 0;           % Left boundary position
G     = h0^2*rhow*g*theta*t0/(12*mu*L0);           % Dimensionless gravity number G (according to your model settings)
t_max = 1000;        % Maximum time
nFrames = 200;       % Number of frames
t_vec = linspace(1e-3, t_max, nFrames);  % Time series (avoiding t=0)

% —— Preload data —— 
% Ensure Sloping_Boundary_Layer.mat is in the current path
load('Sloping_Boundary_Layer.mat','X','H');

% —— Prepare animation output —— 
% (1) MP4 video
v = VideoWriter('thin_film_evolution.mp4','MPEG-4');
v.FrameRate = 20;
open(v);

% (2) GIF file
gif_filename = 'thin_film_evolution.gif';
delayTime = 0.05;  % seconds per frame

% —— Plotting and writing frames —— 
figure('Color','w');
for k = 1:nFrames
    t = t_vec(k);
    
    % Calculate composite solution
    [x,h,~] = Composite_Solution(t, Ll, G);
    
    % Plot
    plot(x, h, 'LineWidth', 2);
    ylim([0, 0.25]);
    xlim([Ll, Ll + (3/2^(2/3))*G^(1/3)*t_max^(1/3)]);
    xlabel('x');
    ylabel('h(x,t)');
    title(sprintf('Thin Film Evolution at t = %.1f', t));
    grid on;
    
    drawnow;
    
    % Write to MP4
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Write to GIF
    img = frame2im(frame);
    [A,map] = rgb2ind(img,256);
    if k == 1
        imwrite(A, map, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end

% —— Finalize —— 
close(v);
fprintf('Animation has been saved as thin_film_evolution.mp4 and thin_film_evolution.gif\n');