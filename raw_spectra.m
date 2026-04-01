clear; close all;

%% -------------------- Config --------------------
ROOT = pwd;

cfg.anglesDeg   = [150, 127.5, 105, 82.5, 60];

% where your HDF5 files live (same as your fit code)
cfg.rawFigRoot  = fullfile(ROOT, 'Figures', 'Raw_spectra');

% thesis-style font size
cfg.fontSize = 18;

ensure_dir(cfg.rawFigRoot);

% -------------------- Materials --------------------
KTO = struct( ...
    'name', 'KTO', ...
    'dataFolder', fullfile(ROOT,'hdf5_data','KTO_data'), ...
    'baseNames', {{ ...
        'KTO100_Qdep_529p4eV_aligned', ...
        'KTO110_Qdep_529p4eV_aligned', ...
        'KTO111_Qdep_529p4eV_aligned' ...
    }} ...
);

STO = struct( ...
    'name', 'STO', ...
    'dataFolder', fullfile(ROOT,'hdf5_data','STO_data'), ...
    'baseNames', {{ ...
        'STO_100_Qdep_aligned', ...
        'STO_110_Qdep_aligned', ...
        'STO_111_Qdep_aligned' ...
    }} ...
);

materials = {KTO, STO};

%% -------------------- Save RAW spectra --------------------
for m = 1:numel(materials)
    mat = materials{m};
    matName = mat.name;

    for b = 1:numel(mat.baseNames)
        baseName = mat.baseNames{b};
        tag      = file_tag_from_basename(baseName);

        h5Path = fullfile(mat.dataFolder, [baseName '.hdf5']);
        if ~exist(h5Path,'file')
            error('HDF5 file not found: %s', h5Path);
        end

        % Output folder: Figures/Raw_spectra/<Material>/<Tag>/
        outDir = fullfile(cfg.rawFigRoot, matName, tag);
        ensure_dir(outDir);

        % Load scans
        aux = load_scans_hdf5_as_aux(h5Path, numel(cfg.anglesDeg));

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

            % Plot raw (blue) with larger fonts
            fig = figure('Visible','off','Color','w');
            ax = axes(fig);
            plot(ax, E, spec, '-', 'LineWidth', 1.2, 'Color', [0 0 1]);
            xlim(ax, [-0.4 1]);   % <-- limit x-axis
            xlabel(ax, 'Energy Loss (eV)', 'FontSize', cfg.fontSize);
            ylabel(ax, 'Intensity (a.u.)', 'FontSize', cfg.fontSize);

            ax.FontSize = cfg.fontSize;
            ax.LineWidth = 1.0;
            ax.Box = 'on';

            % Save
            f = fullfile(outDir, sprintf('raw_%s_%sdeg.png', tag, angle_tag(twoTheta)));
            exportgraphics(fig, f, 'Resolution', 300);
            close(fig);

        end
    end
end

disp('Done: raw spectra saved under Figures/Raw_spectra/<Material>/<Tag>/');

%% ===================== Local functions =====================

function ensure_dir(p)
if ~exist(p,'dir'), mkdir(p); end
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

function s = angle_tag(x)
s = sprintf('%.1f', x);
s = regexprep(s, '0$', '');
s = regexprep(s, '\.$', '');
s = strrep(s, '.', 'p');
end
