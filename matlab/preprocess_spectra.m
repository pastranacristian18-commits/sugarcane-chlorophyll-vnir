function preprocess_spectra()
% PREPROCESS_SPECTRA  Read leaf reflectance files, filter, resample and export matrices.
%
%   Numerically identical to the original 'filtrar_organizar_datos_soca' script
%   (verified: the exported SG+MSC matrix is reproduced to 2e-13).
%
%   ACTUAL processing order:  interpolation -> low-reflectance filter -> MSC -> Savitzky-Golay
%   (The manuscript text says "SG then MSC": correct the text or the code, they must agree.)
%
%   MSC reference = mean spectrum of ALL retained leaves (all campaigns). This is the
%   "global" version; its leakage effect is negligible (see R/scripts/08_msc_sensitivity.R,
%   change < 0.01 % in RMSE). For a fully leak-free pipeline use the SG-only export and
%   re-estimate MSC inside each training fold in R (R/functions.R: pp_msc).
%
%   Outputs (Excel, one row per leaf, columns Dia, Unidad, then 500 wavelengths):
%     DATOS_SOCA_MATRIZ_CRUDA.xlsx   interpolated raw reflectance
%     DATOS_SOCA_SG_SIN_MSC.xlsx     Savitzky-Golay only  (recommended input for fold-wise MSC)
%     DATOS_SOCA_MATRIZ_FINAL.xlsx   MSC (global) followed by Savitzky-Golay

    %% ---- configuration ------------------------------------------------------
    baseDir   = fullfile('data', 'raw_spectra');          % folders 'CAMPAÑA FOLIAR 1' ... 'CAMPAÑA FOLIAR 6'
    outDir    = fullfile('data', 'matlab_export');
    nCamp     = 6;
    wlMin     = 450.064;  wlMax = 820.378;  nBands = 500;
    marker    = '>>>>>Begin Spectral Data<<<<<';
    minRefl   = 1;                                        % discard a spectrum if ANY band < 1 %
    sgWindow  = 11;  sgDegree = 2;
    % files known to be damaged: {campaign, file-name fragment}
    blacklist = {4, 'p4s4vis_reflection'};                % explains the missing Day-4 leaf P4S4
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    wlTarget = linspace(wlMin, wlMax, nBands);

    %% ---- 1. read, filter, resample -------------------------------------------
    keep = struct();
    for c = 1:nCamp
        campDir = fullfile(baseDir, sprintf('CAMPAÑA FOLIAR %d', c));
        if ~exist(campDir, 'dir'), continue; end
        files = dir(fullfile(campDir, '*.txt'));
        for i = 1:numel(files)
            fn = files(i).name;  fl = lower(fn);
            if any(cellfun(@(b) isequal(b{1}, c) && contains(fl, b{2}), num2cell(blacklist, 2))), continue; end
            if ~contains(fl, 'vis') && ~contains(fl, 'vi_'), continue; end
            if contains(fl, 'nir'), continue; end
            tp = regexp(fn, '[pP](\d+)', 'tokens', 'once');  ts = regexp(fn, '[sS](\d+)', 'tokens', 'once');
            if isempty(tp), continue; end
            p = str2double(tp{1});  if isempty(ts), s = 1; else, s = str2double(ts{1}); end
            key = sprintf('C%d_P%d_S%d', c, p, s);  if isfield(keep, key), continue; end
            txt = fileread(fullfile(campDir, fn));  k = strfind(txt, marker);  if isempty(k), continue; end
            data = textscan(strrep(txt(k + numel(marker):end), ',', '.'), '%f %f');
            w = data{1};  r = data{2};  m = w >= wlMin & w <= wlMax & ~isnan(r);  if ~any(m), continue; end
            ri = interp1(w(m), r(m), wlTarget, 'pchip', 'extrap');
            if any(ri < minRefl), continue; end
            keep.(key) = struct('camp', c, 'plot', p, 'leaf', s, 'raw', ri);
        end
    end
    keys = fieldnames(keep);  n = numel(keys);
    if n == 0, error('No valid spectra found.'); end
    Xraw = zeros(n, nBands);  meta = zeros(n, 3);
    for i = 1:n, Xraw(i, :) = keep.(keys{i}).raw;  meta(i, :) = [keep.(keys{i}).camp keep.(keys{i}).plot keep.(keys{i}).leaf]; end

    %% ---- 2. preprocessing -----------------------------------------------------
    ref = mean(Xraw, 1);  Xmsc = Xraw;
    for i = 1:n
        b = polyfit(ref, Xraw(i, :), 1);                  % x = b(2) + b(1) * ref
        if b(1) > 0.1, Xmsc(i, :) = (Xraw(i, :) - b(2)) / b(1); end
    end
    sg = @(M) cell2mat(arrayfun(@(i) smoothdata(M(i, :), 'sgolay', sgWindow, 'Degree', sgDegree), (1:size(M, 1))', 'UniformOutput', false));
    Xfinal = sg(Xmsc);   XsgOnly = sg(Xraw);

    %% ---- 3. export ----------------------------------------------------------------
    [~, ord] = sortrows(meta);
    T = table(arrayfun(@(c) sprintf('Dia %d', c), meta(ord, 1), 'UniformOutput', false), ...
              arrayfun(@(p, s) sprintf('P%dS%d', p, s), meta(ord, 2), meta(ord, 3), 'UniformOutput', false), 'VariableNames', {'Dia', 'Unidad'});
    names = cellstr(string(round(wlTarget, 3)));
    exportMat = @(M, f) writetable([T, array2table(M(ord, :), 'VariableNames', names)], fullfile(outDir, f));
    exportMat(Xraw,   'DATOS_SOCA_MATRIZ_CRUDA.xlsx');
    exportMat(XsgOnly, 'DATOS_SOCA_SG_SIN_MSC.xlsx');
    exportMat(Xfinal, 'DATOS_SOCA_MATRIZ_FINAL.xlsx');
    fprintf('Done: %d leaves exported to %s\n', n, outDir);
end
