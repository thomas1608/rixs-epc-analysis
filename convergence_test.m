clear; close all;

%% -------------------- Inputs --------------------
g1 = 30; g2 = 20; g3 = 10;

mFix  = 40;
tolRel = 1e-2;
W      = 5;

KTO_phononE = [0.018; 0.051; 0.107];
STO_phononE = [0.022; 0.060; 0.100];

channels = [
    1 0 0; 0 1 0; 0 0 1;
    2 0 0; 0 2 0; 0 0 2; 1 1 0; 1 0 1; 0 1 1;
    3 0 0; 0 3 0; 0 0 3; 1 1 1; 2 1 0; 2 0 1; 0 2 1; 1 2 0; 1 0 2; 0 1 2
];

%% -------------------- Output folder --------------------
ROOT = pwd;
outRoot = fullfile(ROOT, 'Figures', 'Convergence_test');
ensure_dir(outRoot);

%% -------------------- Run --------------------
mMax = choose_mmax([g1 g2 g3], mFix);

run_material('KTO', KTO_phononE, channels, g1,g2,g3, mFix, mMax, tolRel, W, outRoot);
run_material('STO', STO_phononE, channels, g1,g2,g3, mFix, mMax, tolRel, W, outRoot);

disp(['Done. Figures saved under: ' outRoot]);

%% ===================== Local functions =====================

function ensure_dir(p)
if ~exist(p,'dir'), mkdir(p); end
end

function mMax = choose_mmax(g, mFix)
gmax = max(g);
mMax = max(mFix + 10, ceil(gmax + 6*sqrt(max(gmax,1)) + 25));
mMax = min(mMax, 300);
end

function run_material(matName, omegas, channels, g1,g2,g3, mFix, mMax, tolRel, W, outRoot)
w1 = omegas(1); w2 = omegas(2); w3 = omegas(3);

matOut = fullfile(outRoot, matName);
ensure_dir(matOut);

run_sweep(matName, 'm1', channels, g1,g2,g3, w1,w2,w3, mFix,mFix,mFix, mMax, tolRel, W, matOut);
run_sweep(matName, 'm2', channels, g1,g2,g3, w1,w2,w3, mFix,mFix,mFix, mMax, tolRel, W, matOut);
run_sweep(matName, 'm3', channels, g1,g2,g3, w1,w2,w3, mFix,mFix,mFix, mMax, tolRel, W, matOut);
end

function run_sweep(matName, whichM, channels, g1,g2,g3, w1,w2,w3, m1Fix,m2Fix,m3Fix, mMax, tolRel, W, outDir)
m_scan = 0:mMax;
nCh = size(channels,1);

Ich  = zeros(numel(m_scan), nCh);
Itot = zeros(numel(m_scan), 1);

for k = 1:numel(m_scan)
    mv = m_scan(k);
    switch whichM
        case 'm1', m1 = mv;    m2 = m2Fix; m3 = m3Fix;
        case 'm2', m1 = m1Fix; m2 = mv;    m3 = m3Fix;
        case 'm3', m1 = m1Fix; m2 = m2Fix; m3 = mv;
        otherwise, error('Unknown whichM: %s', whichM);
    end

    s = 0;
    for c = 1:nCh
        n1 = channels(c,1); n2 = channels(c,2); n3 = channels(c,3);
        val = I(n1,n2,n3, m1,m2,m3, w1,w2,w3, g1,g2,g3);
        Ich(k,c) = val;
        s = s + val;
    end
    Itot(k) = s;
end

[m_conv, didConverge] = find_plateau_m(Itot, m_scan, tolRel, W);

% ----- Plot -----
fig = figure('Color','w'); hold on; grid on;

for c = 1:nCh
    plot(m_scan, Ich(:,c), 'LineWidth', 1);
end
plot(m_scan, Itot, 'k-', 'LineWidth', 2);

if didConverge
    xline(m_conv, 'r--', sprintf('m_{conv}=%d', m_conv), ...
        'LabelVerticalAlignment','bottom', 'LineWidth', 1.5);
else
    xline(m_scan(end), 'r--', sprintf('no conv up to %d', m_scan(end)), ...
        'LabelVerticalAlignment','bottom', 'LineWidth', 1.5);
end

switch whichM
    case 'm1', xlab = 'm_1'; fixedStr = sprintf('(m_2,m_3)=(%d,%d)', m2Fix, m3Fix);
    case 'm2', xlab = 'm_2'; fixedStr = sprintf('(m_1,m_3)=(%d,%d)', m1Fix, m3Fix);
    case 'm3', xlab = 'm_3'; fixedStr = sprintf('(m_1,m_2)=(%d,%d)', m1Fix, m2Fix);
end

xlabel(sprintf('%s  %s', xlab, fixedStr));
ylabel('Intensity |A|^2');

% ----- Font sizes -----
ax = gca;
ax.FontSize = 18;          % tick labels
ax.LineWidth = 1.2;

xlabel(ax.XLabel.String, 'FontSize', 18);
ylabel(ax.YLabel.String, 'FontSize', 18);

% ----- Save -----
base = sprintf('%s_%s_convergence_mFix%d_tol%.0e_W%d', matName, whichM, m1Fix, tolRel, W);
pngPath = fullfile(outDir, [base '.png']);
figPath = fullfile(outDir, [base '.fig']);

exportgraphics(fig, pngPath, 'Resolution', 200);
savefig(fig, figPath);
close(fig);

fprintf('[SAVE] %s\n', pngPath);
end

function [mSat, didConverge] = find_plateau_m(Ivec, m_scan, tolRel, W)
Ivec = Ivec(:);
N = numel(Ivec);
mSat = m_scan(end);
didConverge = false;
if N < (W+1), return; end

eps0 = 1e-12;

for i = 1:(N-W)
    a = Ivec(i);
    b = Ivec(i+W);

    if abs(b) < eps0 && abs(a) < eps0
        continue
    end

    rel = abs(b - a) / max(abs(b), eps0);
    if rel < tolRel
        mSat = m_scan(i+W);
        didConverge = true;
        return;
    end
end
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

A = 0;
for m1 = 0:m1_values
    base1 = z + w1*(g1-m1);
    c1 = C1(m1+1);
    for m2 = 0:m2_values
        base12 = base1 + w2*(g2-m2) + w3*g3;
        A = A + c1 * C2(m2+1) * sum(C3 ./ (base12 - w3*m3_idx));
    end
end

intensity = abs(A)^2;
end

function B = B_nm_g(n,m,g)
persistent FACT
if isempty(FACT)
    FACTMAX = 400;
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
