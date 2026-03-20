% 10 October 2025 experimental nevis and ice velocity run
% this example just tests the velocity calculation, no hydrology

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis/src';   % code directory
addpath(oo.code);       % add path to code

%% parameters
[pd,oo] = nevis_defaults([],oo);
% [ put non-default parameters and options here ]
[ps,pp] = nevis_nondimension(pd,[],oo);

%% grid and geometry
x = linspace(0,(50000/ps.x),101); y = linspace(0,(10000/ps.x),21);
gg = nevis_grid(x,y,oo);
b = (0/ps.z)*gg.nx;
s = (1500/ps.z)*max(1-gg.nx/max(max(gg.nx)),0).^(1/2);

%% mask grid
% gg = nevis_mask(gg,[]); 
gg = nevis_mask(gg,find(s-b<=0)); 
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

%% plot grid
% nevis_plot_grid(gg);

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden 
vv.hs = (0.5/ps.hs)*ones(gg.nIJ,1);         % 50cm thick sheet

%% add parameters and boundary labels for ice velocity
pd.n_glen = 1;
eps = 0.01; pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation;
% pd.n_glen = 3;
% pd.A_glen = pd.A;
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity

gg = nevis_label_ice(gg,oo); % add boundary labels needed for ice velocity 
figure(1); clf; nevis_plot_grid_ice(gg); 

%% initial ice velocity
N = ones(gg.nIJ,1); 
u = zeros(gg.eIJ,1); 
v = zeros(gg.fIJ,1); 
u(gg.eout2) = NaN; 
v(gg.fout2) = NaN;
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
[tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
figure(2); clf; imagesc(gg.ex(:,1),gg.ey(1,:),reshape(u,gg.eI,gg.eJ)'); colorbar;

%% strip of zero effective pressure near middle of domain    
N = ones(gg.nIJ,1); N(gg.nx>2 & gg.nx<3) = 0; 
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
[tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
figure(3); clf; 
    plot((ps.x/1e3)*gg.ex(:,1),(ps.u*pd.ty)*u(1:gg.eI));
    xlabel('x [km]'); ylabel('u [m/y]');
figure(4); clf; 
    plot((ps.x/1e3)*gg.nx(:,1),(ps.eta*ps.u/ps.x/1e6)*Txx(1:gg.nI),'-');
    xlabel('x [km]'); ylabel('\tau_{xx} [MPa]');

%% patch of zero effective pressure near middle of domain
N = ones(gg.nIJ,1); N(((gg.nx-2.5).^2+(gg.ny-0.5).^2).^(1/2)<0.3) = 0; 
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
[tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
figure(3); clf; 
    U = ((gg.nmeanx(:,gg.es2)*u(gg.es2)).^2+(gg.nmeany(:,gg.fs2)*v(gg.fs2)).^2).^(1/2); U(gg.nout) = NaN;
    imagesc((ps.x/1e3)*gg.ex(:,1),(ps.x/1e3)*gg.ey(1,:),(ps.u*pd.ty)*reshape(U,gg.nI,gg.nJ)'); colorbar;
    xlabel('x [km]'); ylabel('y [km]'); 
figure(4); clf; 
    Txx(gg.nout) = NaN; % don't plot meaningless values outside of domain
    imagesc((ps.x/1e3)*gg.nx(:,1),(ps.x/1e3)*gg.ny(1,:),(ps.eta*ps.u/ps.x/1e6)*reshape(Txx,gg.nI,gg.nJ)'); colorbar;
    xlabel('x [km]'); ylabel('y [km]'); 
  

