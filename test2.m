alphah = pd.k_s*ps.h^(pd.alpha_s-1)*ps.phi^pd.beta_s/ps.hb/ps.x^(pd.beta_s+1);
alphas = pd.k_c*ps.S^(pd.alpha_c-1)*ps.phi^(pd.beta_c)/ps.hb/ps.x^(pd.beta_c+2);
disp(['alpha_s = ' num2str(alphas) ' per day, alpha_h = ' num2str(alphah)]);
% disp(['alpha_s = ' num2str(alphas*10*pd.td) ' per day, alpha_h = ' num2str(alphah*ps.h*pd.td) ' per day']);

% alphah = pd.k_s*ps.phi^pd.beta_s/ps.hb/ps.x^(pd.beta_s+1);
% alphas = pd.k_c*ps.phi^(pd.beta_c)/ps.hb/ps.x^(pd.beta_c+2);
% disp(['alpha_s = ' num2str(alphas*(10)^(pd.alpha_c)*pd.td) ' per day, alpha_h = ' num2str(alphah*ps.h^(pd.alpha_s)*pd.td) ' per day']);
alphah = 1/(10*pd.td)/ps.h;
alphas = 1/(10*pd.td)/ps.S;
disp(['alpha_s = ' num2str(alphas) ' per day, alpha_h = ' num2str(alphah)]);