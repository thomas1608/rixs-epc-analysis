
clear; close all;

%% -------------------- Config --------------------
ROOT = pwd;

cfg.fitWin      = [-0.1, 0.25];               % eV
cfg.anglesDeg   = [150, 127.5, 105, 82.5, 60];

cfg.figRoot             = fullfile(ROOT, 'Figures');
cfg.figRootElasticSub   = fullfile(ROOT, 'Figures', 'Elastic_subtracted');

ensure_dir(cfg.figRoot);
ensure_dir(cfg.figRootElasticSub);

% -------------------- Materials --------------------
KTO = struct( ...
    'dataFolder', fullfile(ROOT,'hdf5_data','KTO_data'), ...
    'baseNames', {{ ...
        'KTO100_Qdep_529p4eV_aligned', ...
        'KTO110_Qdep_529p4eV_aligned',...
        'KTO111_Qdep_529p4eV_aligned' ...
    }}, ...
    'phononE', [0.018; 0.051; 0.107], ...
    'mu', 0.475, ...
    'fwhm', 0.023, ...
    'fwhm_STO100', [], ...
    'p0', [30, 20, 10, 6, 0, 0.05], ...
    'lb', [0,  0,  0, 0, -0.02, 0], ...
    'ub', [60, 60, 30, inf, 0.02, inf], ...
    'aA', 4.03, ...
    'Einc', 529.4 ...
);

STO = struct( ...
    'dataFolder', fullfile(ROOT,'hdf5_data','STO_data'), ...
    'baseNames', {{ ...
        'STO_100_Qdep_aligned', ...
        'STO_110_Qdep_aligned', ...
        'STO_111_Qdep_aligned' ...
    }}, ...
    'phononE', [0.022; 0.060; 0.100], ...
    'mu', 0.30, ...
    'fwhm', 0.027, ...
    'fwhm_STO100', 0.025, ...
    'p0', [30, 20, 10, 6, 0, 0.05], ...
    'lb', [0,  0,  0, 0, -0.02, 0], ...
    'ub', [60, 60, 30, inf, 0.02, inf], ...
    'aA', 3.903, ...
    'Einc', 531.0 ...
);

materials = {struct('name','KTO','cfg',KTO), struct('name','STO','cfg',STO)};

%% -------------------- Results prealloc --------------------
varNames = { ...
    'Material','Tag','BaseName','H5File','AngleDeg','RLU', ...
    'g1','g1_plus','g1_minus','g2','g2_plus','g2_minus','g3','g3_plus','g3_minus', ...
    'M1','M1_plus','M1_minus','M2','M2_plus','M2_minus','M3','M3_plus','M3_minus', ...
    'resnorm','par1','par2','par3','par4','par5','par6','res_fwhm_used' ...
};

nRowsExpected = 0;
for m = 1:numel(materials)
    nRowsExpected = nRowsExpected + numel(materials{m}.cfg.baseNames) * numel(cfg.anglesDeg);
end

rows = cell(nRowsExpected, numel(varNames));
rowi = 0;

cfg.timing.doPrintPerAngle = true;
timingVarNames = {'Material','Tag','AngleDeg','t_fit_s','resnorm'};
TimingRows = cell(nRowsExpected, numel(timingVarNames));
timingi = 0;

t_total_all = tic;

%% -------------------- Run fits --------------------
for m = 1:numel(materials)
    matName = materials{m}.name;
    MC      = materials{m}.cfg;

    t_material = tic;
    nFitsThisMat = 0;
    tFitsThisMat = [];

    for b = 1:numel(MC.baseNames)
        baseName = MC.baseNames{b};
        tag      = file_tag_from_basename(baseName);
        hkl      = hkl_from_tag(tag);

        t_tag = tic;
        tFitsThisTag = [];

        h5Path = fullfile(MC.dataFolder, [baseName '.hdf5']);
        if ~exist(h5Path,'file'), error('HDF5 file not found: %s', h5Path); end

        figDir    = fullfile(cfg.figRoot, matName, tag);
        figDirNoE = fullfile(cfg.figRootElasticSub, matName, tag);
        ensure_dir(figDir);
        ensure_dir(figDirNoE);

        aux = load_scans_hdf5_as_aux(h5Path, numel(cfg.anglesDeg));
        p0  = MC.p0;

        for ii = 1:numel(cfg.anglesDeg)
            twoTheta = cfg.anglesDeg(ii);
            E    = aux(:, 2*ii-1);
            spec = aux(:, 2*ii);

            isSTO100 = strcmp(matName,'STO') && strcmp(tag,'STO100');
            if ~isSTO100
                spec = spec / 1000;
            else
                spec = subtract_bg_sto100(E, spec, -0.20, -0.10);
            end

            iMin = find(E > cfg.fitWin(1), 1, 'first'); if isempty(iMin), iMin = 1; end
            iMax = find(E > cfg.fitWin(2), 1, 'first'); if isempty(iMax), iMax = numel(E); end

            rlu = two_theta_to_rlu(twoTheta, MC.aA, MC.Einc);

            w_use = MC.fwhm;
            if isSTO100 && ~isempty(MC.fwhm_STO100)
                w_use = MC.fwhm_STO100;
            end

            t_angle = tic;
            [~, ~, pFit, ci, resnorm] = completeFit_threePhonons( ...
                E(iMin:iMax), spec(iMin:iMax), p0, ...
                MC.lb, MC.ub, w_use, MC.mu, MC.phononE, ...
                figDir, figDirNoE, tag, matName, twoTheta, rlu);
            t_fit = toc(t_angle);

            nFitsThisMat = nFitsThisMat + 1;
            tFitsThisMat(end+1,1) = t_fit;
            tFitsThisTag(end+1,1) = t_fit;

            timingi = timingi + 1;
            TimingRows(timingi,:) = {matName, tag, twoTheta, t_fit, resnorm};

            if cfg.timing.doPrintPerAngle
                fprintf('[TIMING] %s %-6s 2θ=%5.1f° : %.2f s (resnorm=%.3g)\n', ...
                    matName, tag, twoTheta, t_fit, resnorm);
            end

            % g and propagated M (meV)
            g     = pFit(1:3).';
            ci_g  = ci(1:3,:);
            g_min = g - ci_g(:,1);
            g_pls = ci_g(:,2) - g;

            omegas = MC.phononE;
            M_meV  = 1000 * (omegas .* sqrt(max(g,0)));

            eps_ = 1e-300;
            coef = omegas ./ (2*sqrt(max(g,eps_)));
            M_min = 1000 * coef .* g_min;
            M_pls = 1000 * coef .* g_pls;

            rowi = rowi + 1;
            rows(rowi,:) = { ...
                matName, tag, baseName, h5Path, twoTheta, rlu, ...
                g(1), g_pls(1), g_min(1), ...
                g(2), g_pls(2), g_min(2), ...
                g(3), g_pls(3), g_min(3), ...
                M_meV(1), M_pls(1), M_min(1), ...
                M_meV(2), M_pls(2), M_min(2), ...
                M_meV(3), M_pls(3), M_min(3), ...
                resnorm, ...
                pFit(1), pFit(2), pFit(3), pFit(4), pFit(5), pFit(6), ...
                w_use ...
            };

        end

        fprintf('[TIMING] %s %-6s (%d angles): total %.2f s | mean %.2f s | median %.2f s | max %.2f s\n', ...
            matName, tag, numel(cfg.anglesDeg), toc(t_tag), ...
            mean(tFitsThisTag), median(tFitsThisTag), max(tFitsThisTag));
    end

    fprintf('\n[TIMING] %s (all tags): total %.2f s | fits=%d | mean %.2f s | median %.2f s | max %.2f s\n\n', ...
        matName, toc(t_material), nFitsThisMat, mean(tFitsThisMat), median(tFitsThisMat), max(tFitsThisMat));
end

%% -------------------- Save tables --------------------
rows    = rows(1:rowi,:);
Results = cell2table(rows, 'VariableNames', varNames);


Compact = Results(:, { ...
    'Material','Tag','AngleDeg','RLU', ...
    'g1','g1_plus','g1_minus','g2','g2_plus','g2_minus','g3','g3_plus','g3_minus', ...
    'M1','M1_plus','M1_minus','M2','M2_plus','M2_minus','M3','M3_plus','M3_minus', ...
    'resnorm','res_fwhm_used' ...
});

%% -------------------- Print summaries --------------------
disp(' ');
disp('==================== Qdep Fit Summary (g and M with 95% CI) ====================');

Results = sortrows(Results, {'Material','Tag','AngleDeg'}, {'ascend','ascend','descend'});

%% -------------------- Plot M vs scattering angle (save figures) --------------------
outDirM = fullfile(cfg.figRoot, 'Summary_M_vs_angle');
ensure_dir(outDirM);

plot_M_vs_angle_alltags(Results, outDirM);
disp(['M-vs-angle figures saved under:     ' outDirM]);

disp(' ');
disp('--- Coupling constants g (95% CI) ---');
pretty_print_g_table(Results(:, {'Material','Tag','AngleDeg','RLU','g1','g1_plus','g1_minus','g2','g2_plus','g2_minus','g3','g3_plus','g3_minus','resnorm','res_fwhm_used'}));

disp(' ');
disp('--- EPC M (meV) (propagated from g CIs) ---');
pretty_print_M_table(Results(:, {'Material','Tag','AngleDeg','RLU','M1','M1_plus','M1_minus','M2','M2_plus','M2_minus','M3','M3_plus','M3_minus','res_fwhm_used'}));

disp(' ');
disp(['Figures saved under:               ' cfg.figRoot]);
disp(['Elastic-subtracted figures under:  ' cfg.figRootElasticSub]);

%% -------------------- Timing summary --------------------
TimingRows = TimingRows(1:timingi,:);
Timing     = cell2table(TimingRows, 'VariableNames', timingVarNames);

t_total = toc(t_total_all);

disp(' ');
disp('==================== TIMING SUMMARY ====================');
fprintf('Total runtime (both materials): %.2f s (%.2f min)\n', t_total, t_total/60);

disp(' ');
disp('--- Per material (per-angle fit time) ---');
disp(groupsummary(Timing, "Material", ["mean","median","max","sum"], "t_fit_s"));

disp(' ');
disp('--- Per tag (per-angle fit time) ---');
disp(groupsummary(Timing, ["Material","Tag"], ["mean","median","max","sum"], "t_fit_s"));

disp('Done.');

%% ===================== Local functions =====================

function ensure_dir(p)
if ~exist(p,'dir'), mkdir(p); end
end

function pretty_print_g_table(T)
fprintf('%-4s %-6s %7s %9s %7s | %-22s %-22s %-22s | %9s\n', ...
    'Mat','Tag','2θ','r.l.u.','FWHM','g1 (+/-)','g2 (+/-)','g3 (+/-)','resnorm');
fprintf('%s\n', repmat('-',1,4+1+6+1+7+1+9+1+7+3+22*3+3+9));
for i = 1:height(T)
    fprintf('%-4s %-6s %7.1f %9.4f %7.3f | ', ...
        string(T.Material(i)), string(T.Tag(i)), T.AngleDeg(i), T.RLU(i), T.res_fwhm_used(i));
    fprintf('%.3f (+%.3f/-%.3f) ', T.g1(i), T.g1_plus(i), T.g1_minus(i));
    fprintf('%.3f (+%.3f/-%.3f) ', T.g2(i), T.g2_plus(i), T.g2_minus(i));
    fprintf('%.3f (+%.3f/-%.3f) ', T.g3(i), T.g3_plus(i), T.g3_minus(i));
    fprintf('| %9.3g\n', T.resnorm(i));
end
end

function pretty_print_M_table(T)
fprintf('%-4s %-6s %7s %9s %7s | %-22s %-22s %-22s\n', ...
    'Mat','Tag','2θ','r.l.u.','FWHM','M1 (+/-) meV','M2 (+/-) meV','M3 (+/-) meV');
fprintf('%s\n', repmat('-',1,4+1+6+1+7+1+9+1+7+3+22*3));
for i = 1:height(T)
    fprintf('%-4s %-6s %7.1f %9.4f %7.3f | ', ...
        string(T.Material(i)), string(T.Tag(i)), T.AngleDeg(i), T.RLU(i), T.res_fwhm_used(i));
    fprintf('%.2f (+%.2f/-%.2f) ', T.M1(i), T.M1_plus(i), T.M1_minus(i));
    fprintf('%.2f (+%.2f/-%.2f) ', T.M2(i), T.M2_plus(i), T.M2_minus(i));
    fprintf('%.2f (+%.2f/-%.2f)\n', T.M3(i), T.M3_plus(i), T.M3_minus(i));
end
end

function aux = load_scans_hdf5_as_aux(h5Path, nScans)
info  = h5info(h5Path);
names = {info.Datasets.Name};
names = names(cellfun(@(s) ~isempty(regexp(s,'^\d+$','once')), names));
if isempty(names), error('No numeric scan datasets in %s', h5Path); end

idx = cellfun(@str2double, names);
[~, ord] = sort(idx);
names = names(ord);

if numel(names) < nScans
    error('Need %d scans but file has %d', nScans, numel(names));
end

d0  = ensure_N_by_2(h5read(h5Path, ['/' names{1}]), names{1});
N   = size(d0,1);
aux = zeros(N, 2*nScans);
aux(:,1:2) = d0;

for i = 2:nScans
    d = ensure_N_by_2(h5read(h5Path, ['/' names{i}]), names{i});
    if size(d,1) ~= N, error('Dataset %s row count mismatch', names{i}); end
    aux(:, 2*i-1:2*i) = d;
end
end

function d = ensure_N_by_2(d, name)
if ~ismatrix(d), error('Dataset %s not 2D', name); end
if size(d,2) == 2
    return
elseif size(d,1) == 2
    d = d.';
else
    error('Dataset %s is not (N×2) or (2×N)', name);
end
end

function y = subtract_bg_sto100(energy, y, lo, hi)
m = (energy >= lo) & (energy <= hi);
if any(m)
    y = y - mean(y(m));
end
end

function tag = file_tag_from_basename(baseName)
parts = strsplit(baseName, '_');
if numel(parts) >= 2 && all(isstrprop(parts{2},'digit'))
    tag = [parts{1} parts{2}];
else
    tag = parts{1};
end
end

function hkl = hkl_from_tag(tag)
d = regexp(tag, '\d', 'match');
d = [d{:}];
if numel(d) ~= 3, error('Could not extract hkl from tag: %s', tag); end
hkl = [str2double(d(1)), str2double(d(2)), str2double(d(3))];
end

function s = angle_tag(x)
s = sprintf('%.1f', x);
s = regexprep(s, '0$', '');
s = regexprep(s, '\.$', '');
s = strrep(s, '.', 'p');
end

function rlu = two_theta_to_rlu(twoThetaDeg, latticeA_Ang, incidentE_eV)
theta   = deg2rad(twoThetaDeg / 2.0);
lambdaA = 12398.4 / incidentE_eV;
Q       = (4.0*pi/lambdaA) * sin(theta);   % Å^-1
rlu     = (latticeA_Ang/(2.0*pi)) * Q;     % Q / (2π/a)
end


function [fL, fR] = get_elastic_fwhm(matName, tag, twoTheta, defaultL, defaultR)
fL = defaultL; fR = defaultR;
if strcmp(matName,'STO') && abs(twoTheta - 150.0) < 1e-12 && any(strcmp(tag, {'STO100','STO110','STO111'}))
    fL = defaultL + 0.003;
    fR = defaultR - 0.001;
end
end

function [fitting, data_no_elastic, par, ci, resnorm] = completeFit_threePhonons( ...
    energy, scan, par0, lb, ub, w, mu, omegas, saveDir, saveDirNoEl, fileTag, matName, twoTheta, rlu)

w1 = omegas(1); w2 = omegas(2); w3 = omegas(3);
ms1 = 50; ms2 = 40; ms3 = 25;

sigma = w/(2*sqrt(2*log(2)));
norm_factor = I(1,0,0,ms1,ms2,ms3,w1,w2,w3,15,15,15);

pv_ph = @(amp,xc,x) pseudoVoigtAsymmetric([amp, xc, w+0.001, w-0.001, mu], x);

% elastic can have per-dataset override
defaultL = w + 0.001;
defaultR = w - 0.001;
[fwhmL_el, fwhmR_el] = get_elastic_fwhm(matName, fileTag, twoTheta, defaultL, defaultR);
pv_el = @(amp,xc,x) pseudoVoigtAsymmetric([amp, xc, fwhmL_el, fwhmR_el, mu], x);

v = @(amp,xc,x) amp/2 .* exp(-0.5*((x-xc)./(sigma*2)).^2);

modelFun = @(k,x) k(6)/norm_factor*( ...
    I(1,0,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*pv_ph(1, w1+k(5), x) + ...
    I(0,1,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*pv_ph(1, w2+k(5), x) + ...
    I(0,0,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*pv_ph(1, w3+k(5), x) + ...
    I(2,0,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 2*w1+k(5), x) + ...
    I(0,2,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 2*w2+k(5), x) + ...
    I(0,0,2,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 2*w3+k(5), x) + ...
    I(1,1,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w1+w2+k(5), x) + ...
    I(1,0,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w1+w3+k(5), x) + ...
    I(0,1,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w2+w3+k(5), x) + ...
    I(3,0,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 3*w1+k(5), x) + ...
    I(0,3,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 3*w2+k(5), x) + ...
    I(0,0,3,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 3*w3+k(5), x) + ...
    I(1,1,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w1+w2+w3+k(5), x) + ...
    I(2,1,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 2*w1+w2+k(5), x) + ...
    I(2,0,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 2*w1+w3+k(5), x) + ...
    I(0,2,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, 2*w2+w3+k(5), x) + ...
    I(1,2,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w1+2*w2+k(5), x) + ...
    I(1,0,2,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w1+2*w3+k(5), x) + ...
    I(0,1,2,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v(1, w2+2*w3+k(5), x) ) + ...
    pv_el(k(4), k(5), x);

opts = optimoptions('lsqcurvefit', ...
    'Display','off', 'MaxFunctionEvaluations', 5000, 'MaxIterations', 2000);

exitflag = 0;
while exitflag == 0
    [par,resnorm,residual,exitflag,~,~,jacobian] = lsqcurvefit(modelFun, par0, energy, scan, lb, ub, opts);
    par0 = par;
end

fitting = modelFun(par,energy);
ci      = nlparci(par,residual,'jacobian',jacobian);

k4 = par(4); k5 = par(5);
elastic_fit = pv_el(k4, k5, energy);
data_no_elastic  = scan - elastic_fit;
model_no_elastic = fitting - elastic_fit;

% --- fixed component colors (same for both figures) ---
C.fit   = [0.85 0.00 0.00];   % red
C.el    = [0.93 0.69 0.13];   % yellow-ish
C.ph1   = [0.89 0.34 0.09];   % orange
C.ph2   = [0.60 0.10 1.00];   % purple
C.ph3   = [0.30 0.75 0.93];   % cyan
C.ov2   = [0.35 0.35 0.35];   % dark gray (degree-2 overtones)
C.ov3   = [0.63 0.10 0.22];   % dark red (degree-3 overtones)
C.res   = [0.70 0.70 0.70];   % light gray


% Plot 1: full fit
fig = figure('Visible','off','Color','w');
energy2 = linspace(energy(1),energy(end),5*numel(energy))';
plot(energy-k5, scan, 'ko-','markersize',4,'linewidth',0.7, ...
     'DisplayName',sprintf('2\\theta = %.1f^\\circ (%.4f r.l.u.)', twoTheta, rlu)); hold on

plot(energy-k5, fitting, 'LineWidth',2, 'Color',C.fit, 'DisplayName','Fitting');
plot(energy2-k5, pv_el(k4,k5,energy2), '-', 'LineWidth',2, 'Color',C.el, 'DisplayName','Elastic');

% components (same as before)
k1 = par(1); k2 = par(2); k3 = par(3); k6 = par(6)/norm_factor;
I100 = k6*I(1,0,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I010 = k6*I(0,1,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I001 = k6*I(0,0,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I200 = k6*I(2,0,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I020 = k6*I(0,2,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I002 = k6*I(0,0,2,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I110 = k6*I(1,1,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I101 = k6*I(1,0,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I011 = k6*I(0,1,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I300 = k6*I(3,0,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I030 = k6*I(0,3,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I003 = k6*I(0,0,3,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I111 = k6*I(1,1,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I210 = k6*I(2,1,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I201 = k6*I(2,0,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I021 = k6*I(0,2,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I120 = k6*I(1,2,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I102 = k6*I(1,0,2,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I012 = k6*I(0,1,2,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);

plot(energy2-k5, I100*pv_ph(1,w1+k5,energy2), '-', 'LineWidth',2, 'Color',C.ph1, 'DisplayName','Phonon 1');
plot(energy2-k5, I010*pv_ph(1,w2+k5,energy2), '-', 'LineWidth',2, 'Color',C.ph2, 'DisplayName','Phonon 2');
plot(energy2-k5, I001*pv_ph(1,w3+k5,energy2), '-', 'LineWidth',2, 'Color',C.ph3, 'DisplayName','Phonon 3');

plot(energy2-k5, I200*v(1,2*w1+k5,energy2)+I300*v(1,3*w1+k5,energy2), '--', 'LineWidth',2, 'Color',C.ph1, 'DisplayName','2*Ph1 + 3*Ph1');
plot(energy2-k5, I020*v(1,2*w2+k5,energy2)+I030*v(1,3*w2+k5,energy2), '--', 'LineWidth',2, 'Color',C.ph2, 'DisplayName','2*Ph2 + 3*Ph2');
plot(energy2-k5, I002*v(1,2*w3+k5,energy2)+I003*v(1,3*w3+k5,energy2), '--', 'LineWidth',2, 'Color',C.ph3, 'DisplayName','2*Ph3 + 3*Ph3');

plot(energy2-k5, I110*v(1,w1+w2+k5,energy2)+I101*v(1,w1+w3+k5,energy2)+I011*v(1,w2+w3+k5,energy2), ...
     '--','LineWidth',2,'Color',C.ov2,'DisplayName','other degree-2');

plot(energy2-k5, I111*v(1,w1+w2+w3+k5,energy2)+I210*v(1,2*w1+w2+k5,energy2)+I201*v(1,2*w1+w3+k5,energy2)+ ...
              I021*v(1,2*w2+w3+k5,energy2)+I120*v(1,w1+2*w2+k5,energy2)+I102*v(1,w1+2*w3+k5,energy2)+ ...
              I012*v(1,w2+2*w3+k5,energy2), ...
     '--','LineWidth',2,'Color',C.ov3,'DisplayName','other degree-3');

plot(energy-k5, scan - fitting, 'LineWidth',1.5, 'Color',C.res, 'DisplayName','residual');


xlabel('Energy Loss (eV)'); xlim([-0.1 0.25])
yl = ylim; plot([0,0],[yl(1),yl(2)],'k','linewidth',1,'HandleVisibility','off');
ylabel('Intensity (a.u.)'); legend('Location','northeast');

f1 = fullfile(saveDir, sprintf('fitting_%s_%sdeg', fileTag, angle_tag(twoTheta)));
exportgraphics(fig, [f1 '.png'], 'Resolution', 200);
close(fig);

% Plot 2: elastic-subtracted
fig2 = figure('Visible','off','Color','w');
plot(energy-k5, data_no_elastic, 'ko-','markersize',4,'linewidth',0.7, ...
     'DisplayName',sprintf('2\\theta = %.1f^\\circ (%.4f r.l.u.)', twoTheta, rlu)); hold on
plot(energy-k5, model_no_elastic, 'LineWidth',2, 'Color',C.fit, 'DisplayName','Fit (no elastic)');

plot(energy2-k5, I100*pv_ph(1,w1+k5,energy2), '-', 'LineWidth',2, 'Color',C.ph1, 'DisplayName','Phonon 1');
plot(energy2-k5, I010*pv_ph(1,w2+k5,energy2), '-', 'LineWidth',2, 'Color',C.ph2, 'DisplayName','Phonon 2');
plot(energy2-k5, I001*pv_ph(1,w3+k5,energy2), '-', 'LineWidth',2, 'Color',C.ph3, 'DisplayName','Phonon 3');

plot(energy2-k5, I200*v(1,2*w1+k5,energy2)+I300*v(1,3*w1+k5,energy2), '--', 'LineWidth',2, 'Color',C.ph1, 'DisplayName','2*Ph1 + 3*Ph1');
plot(energy2-k5, I020*v(1,2*w2+k5,energy2)+I030*v(1,3*w2+k5,energy2), '--', 'LineWidth',2, 'Color',C.ph2, 'DisplayName','2*Ph2 + 3*Ph2');
plot(energy2-k5, I002*v(1,2*w3+k5,energy2)+I003*v(1,3*w3+k5,energy2), '--', 'LineWidth',2, 'Color',C.ph3, 'DisplayName','2*Ph3 + 3*Ph3');

plot(energy2-k5, I110*v(1,w1+w2+k5,energy2)+I101*v(1,w1+w3+k5,energy2)+I011*v(1,w2+w3+k5,energy2), ...
     '--','LineWidth',2,'Color',C.ov2,'DisplayName','other degree-2');

plot(energy2-k5, I111*v(1,w1+w2+w3+k5,energy2)+I210*v(1,2*w1+w2+k5,energy2)+I201*v(1,2*w1+w3+k5,energy2)+ ...
              I021*v(1,2*w2+w3+k5,energy2)+I120*v(1,w1+2*w2+k5,energy2)+I102*v(1,w1+2*w3+k5,energy2)+ ...
              I012*v(1,w2+2*w3+k5,energy2), ...
     '--','LineWidth',2,'Color',C.ov3,'DisplayName','other degree-3');

xlabel('Energy Loss (eV)'); xlim([-0.1 0.25]);
yl = ylim; plot([0,0],[yl(1),yl(2)],'k','linewidth',1,'HandleVisibility','off');
ylabel('Intensity (a.u.)'); legend('Location','northeast');

f2 = fullfile(saveDirNoEl, sprintf('fitting_%s_%sdeg_noelastic', fileTag, angle_tag(twoTheta)));
exportgraphics(fig2, [f2 '.png'], 'Resolution', 200);
close(fig2);
end

function y = pseudoVoigtAsymmetric(par,x)
amp=par(1); xc=par(2); fL=par(3); fR=par(4); W=par(5);
sL=fL/(2*sqrt(2*log(2))); sR=fR/(2*sqrt(2*log(2)));
gL=exp(-0.5*((x-xc)./sL).^2); lL=(fL/2).^2./((x-xc).^2+(fL/2).^2);
gR=exp(-0.5*((x-xc)./sR).^2); lR=(fR/2).^2./((x-xc).^2+(fR/2).^2);
y = zeros(size(x));
m = x < xc;
y(m)  = amp.*((1-W).*gL(m) + W.*lL(m));
y(~m) = amp.*((1-W).*gR(~m) + W.*lR(~m));
end



function intensity = I(n1,n2,n3,m1_values,m2_values,m3_values,w1,w2,w3,g1,g2,g3)
z = 1i*0.150;

B1_nm = zeros(m1_values+1,1); B1_m0 = zeros(m1_values+1,1);
for m1 = 0:m1_values
    B1_nm(m1+1) = B_nm_g(max(n1,m1), min(n1,m1), g1);
    B1_m0(m1+1) = B_nm_g(m1, 0, g1);
end
C1 = B1_nm .* B1_m0;

B2_nm = zeros(m2_values+1,1); B2_m0 = zeros(m2_values+1,1);
for m2 = 0:m2_values
    B2_nm(m2+1) = B_nm_g(max(n2,m2), min(n2,m2), g2);
    B2_m0(m2+1) = B_nm_g(m2, 0, g2);
end
C2 = B2_nm .* B2_m0;

B3_nm = zeros(m3_values+1,1); B3_m0 = zeros(m3_values+1,1);
for m3 = 0:m3_values
    B3_nm(m3+1) = B_nm_g(max(n3,m3), min(n3,m3), g3);
    B3_m0(m3+1) = B_nm_g(m3, 0, g3);
end
C3 = B3_nm .* B3_m0;

m3_idx = (0:m3_values).';

Intensity_number = 0;
for m1 = 0:m1_values
    base1 = z + w1*(g1-m1);
    c1 = C1(m1+1);
    for m2 = 0:m2_values
        base12 = base1 + w2*(g2-m2) + w3*g3;
        s3 = sum(C3 ./ (base12 - w3*m3_idx));
        Intensity_number = Intensity_number + c1 * C2(m2+1) * s3;
    end
end

intensity = abs(Intensity_number)^2;
end


function B = B_nm_g(n,m,g)
persistent FACT
if isempty(FACT)
    FACTMAX = 200;
    FACT = factorial(0:FACTMAX);
end

pow = sqrt(g)^(n-m);

l  = 0:m;
d1 = FACT(l+1);
d2 = FACT(m-l+1);
d3 = FACT(n-m+l+1);

s = sum(((-g).^l .* pow) ./ (d2 .* d1 .* d3));

B = (-1)^n * sqrt(exp(-g)*FACT(n+1)*FACT(m+1)) * s;
end

function plot_M_vs_angle_alltags(Results, outDir)
% Plots M1,M2,M3 vs RLU for each (Material,Tag)
% Includes ALL angles (no filtering) and saves PNG.

mats = unique(string(Results.Material));
for im = 1:numel(mats)
    mat = mats(im);
    Tm  = Results(string(Results.Material)==mat, :);

    tags = unique(string(Tm.Tag));
    for it = 1:numel(tags)
        tag = tags(it);
        T   = Tm(string(Tm.Tag)==tag, :);

        if isempty(T)
            continue;
        end

        % --- sort on x (RLU) ---
        T = sortrows(T, "RLU", "ascend");
        x = T.RLU;

        % Asymmetric errors: errorbar(x,y,neg,pos)
        y1 = T.M1; n1 = T.M1_minus; p1 = T.M1_plus;
        y2 = T.M2; n2 = T.M2_minus; p2 = T.M2_plus;
        y3 = T.M3; n3 = T.M3_minus; p3 = T.M3_plus;

        fig = figure('Visible','off','Color','w');
        hold on; box on;

        errorbar(x, y1, n1, p1, 'o-','LineWidth',1.2,'DisplayName','M1');
        errorbar(x, y2, n2, p2, 'o-','LineWidth',1.2,'DisplayName','M2');
        errorbar(x, y3, n3, p3, 'o-','LineWidth',1.2,'DisplayName','M3');

        xlabel('Momentum transfer (r.l.u.)');
        ylabel('e-ph coupling constant M (meV)');

        % x-limits to match plotted range
        xlim([min(x) max(x)]);
        set(gca,'XDir','normal');
        legend('Location','best');

        % ----- Font sizes -----
        ax = gca;
        ax.FontSize  = 18;   % tick numbers
        ax.LineWidth = 1.2;
        
        xlabel('Momentum transfer (r.l.u.)', 'FontSize', 18);
        ylabel('e-ph coupling constant M (meV)', 'FontSize', 18);
        
        lgd = legend('Location','northeast');
        lgd.FontSize = 14;
        lgd.Box = 'off';
        title(lgd, sprintf('%s ', tag), 'FontSize', 14);


        % Save (now truly M vs RLU)
        fbase = fullfile(outDir, sprintf('M_vs_RLU_%s', tag));
        exportgraphics(fig, [fbase '.png'], 'Resolution', 200);

        close(fig);
    end
end
end
