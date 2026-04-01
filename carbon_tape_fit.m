clear; close all; clc;

%% -------------------- Load carbon tape data --------------------
h5Path = fullfile(pwd,'hdf5_data','carbon_tape_KTO_data','carbon_tape_KTO.hdf5');

info  = h5info(h5Path);
names = {info.Datasets.Name};

% take first dataset
data = h5read(h5Path, ['/' names{1}]);

% ensure Nx2
if size(data,2) ~= 2
    data = data.';
end

E    = data(:,1);
spec = data(:,2);

%% -------------------- Fit window around elastic --------------------
fitWin = [-0.05 0.05];   % adjust if needed

idx = (E >= fitWin(1)) & (E <= fitWin(2));
Efit = E(idx);
Yfit = spec(idx);

%% -------------------- Initial guesses --------------------
amp0 = max(Yfit);
xc0  = Efit(Yfit==max(Yfit));
fwhm0 = 0.03;

%% =========================================================
%% 1) Gaussian fit
%% =========================================================

gauss = @(p,x) p(1)*exp(-4*log(2)*((x-p(2))/p(3)).^2);

p0_g = [amp0, xc0(1), fwhm0];
lb_g = [0, -0.02, 0];
ub_g = [Inf, 0.02, 0.2];

p_g = lsqcurvefit(gauss, p0_g, Efit, Yfit, lb_g, ub_g);

%% =========================================================
%% 2) Symmetric pseudo-Voigt
%% =========================================================
pseudoVoigtSym = @(p,x) p(1)*((1-p(4))* ...
    exp(-4*log(2)*((x-p(2))/p(3)).^2) + ...
    p(4)*( (p(3)/2)^2 ./ ((x-p(2)).^2 + (p(3)/2)^2) ));

p0_pv = [amp0, xc0(1), fwhm0, 0.5];
lb_pv = [0, -0.02, 0, 0];
ub_pv = [Inf, 0.02, 0.2, 1];

p_pv = lsqcurvefit(pseudoVoigtSym, p0_pv, Efit, Yfit, lb_pv, ub_pv);

%% =========================================================
%% 3) Asymmetric pseudo-Voigt
%% =========================================================
pseudoVoigtAsym = @(p,x) asymPV(p,x);

p0_apv = [amp0, xc0(1), fwhm0, fwhm0, 0.5];
lb_apv = [0, -0.02, 0, 0, 0];
ub_apv = [Inf, 0.02, 0.2, 0.2, 1];

p_apv = lsqcurvefit(pseudoVoigtAsym, p0_apv, Efit, Yfit, lb_apv, ub_apv);

%% -------------------- Plot results --------------------
figure; hold on; box on;

%% -------------------- Plot results --------------------
fig = figure('Color','w'); 
hold on; box on;

plot(Efit, Yfit, 'ko', 'DisplayName','Data');

plot(Efit, gauss(p_g,Efit), 'r','LineWidth',2,'DisplayName','Gaussian');
plot(Efit, pseudoVoigtSym(p_pv,Efit), 'b','LineWidth',2,'DisplayName','Symmetric PV');
plot(Efit, pseudoVoigtAsym(p_apv,Efit), 'm','LineWidth',2,'DisplayName','Asymmetric PV');

xlabel('Energy Loss (eV)', 'FontSize', 18);
ylabel('Intensity (a.u.)', 'FontSize', 18);

ax = gca;
ax.FontSize = 16;     % tick labels
ax.LineWidth = 1.2;

legend('Location','best','FontSize',14);

%% -------------------- Save figure --------------------
outDir = fullfile(pwd,'Figures','Carbon_tape');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

exportgraphics(fig, fullfile(outDir,'carbon_tape_fit_comparison.png'), ...
               'Resolution',300);

%% -------------------- Print FWHM values --------------------
fprintf('\nGaussian FWHM      = %.4f eV\n', p_g(3));
fprintf('Symmetric PV FWHM  = %.4f eV\n', p_pv(3));
fprintf('Asymmetric PV FWHM L = %.4f eV\n', p_apv(3));
fprintf('Asymmetric PV FWHM R = %.4f eV\n', p_apv(4));

%% -------------------- Compute R^2 values --------------------

% Total sum of squares
SS_tot = sum((Yfit - mean(Yfit)).^2);

% Residual sum of squares
SS_res_g   = sum((Yfit - gauss(p_g,Efit)).^2);
SS_res_pv  = sum((Yfit - pseudoVoigtSym(p_pv,Efit)).^2);
SS_res_apv = sum((Yfit - pseudoVoigtAsym(p_apv,Efit)).^2);

% R^2
R2_g   = 1 - SS_res_g/SS_tot;
R2_pv  = 1 - SS_res_pv/SS_tot;
R2_apv = 1 - SS_res_apv/SS_tot;

fprintf('\nR^2 values:\n');
fprintf('Gaussian        : %.6f\n', R2_g);
fprintf('Symmetric PV    : %.6f\n', R2_pv);
fprintf('Asymmetric PV   : %.6f\n', R2_apv);

%% =========================================================
%% Asymmetric PV function
%% =========================================================
function y = asymPV(par,x)

amp = par(1);
xc  = par(2);
fL  = par(3);
fR  = par(4);
W   = par(5);

sL = fL/(2*sqrt(2*log(2)));
sR = fR/(2*sqrt(2*log(2)));

gL = exp(-0.5*((x-xc)./sL).^2);
lL = (fL/2).^2 ./ ((x-xc).^2 + (fL/2).^2);

gR = exp(-0.5*((x-xc)./sR).^2);
lR = (fR/2).^2 ./ ((x-xc).^2 + (fR/2).^2);

y = zeros(size(x));

left  = x < xc;
right = ~left;

y(left)  = amp*((1-W).*gL(left)  + W.*lL(left));
y(right) = amp*((1-W).*gR(right) + W.*lR(right));
end
