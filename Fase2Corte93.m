clc;
clear;
close all;

%% ============================================================
% FASE 2
% Preprocesamiento comparativo
%
% Archivo: t2-t25.00.nii
% Corte: 93
%
% Esta fase compara filtros y justifica el uso de mediana 3x3
% como apoyo para detección, sin reemplazar la imagen original.
%% ============================================================

targetFile  = 't2-t25.00.nii';
targetCase  = 't2_t25_00';
targetSlice = 93;

rootFolder = uigetdir(pwd, ...
    'Selecciona la carpeta MRI_t2w_nii');

if isequal(rootFolder, 0)
    error('No seleccionaste carpeta.');
end

phase2Folder = fullfile(rootFolder, ...
    'Resultados_t2_t25_00_Corte93', ...
    'FaseII_Preprocesamiento');

if ~exist(phase2Folder, 'dir')
    mkdir(phase2Folder);
end

%% Buscar archivo NIfTI

niiPath = fullfile(rootFolder, targetFile);

if ~exist(niiPath, 'file')
    files = dir(fullfile(rootFolder, '**', targetFile));

    if isempty(files)
        error('No se encontró el archivo %s dentro de la carpeta seleccionada.', targetFile);
    end

    niiPath = fullfile(files(1).folder, files(1).name);
end

fprintf('Archivo encontrado:\n%s\n', niiPath);

%% Leer volumen

info = niftiinfo(niiPath);
volume = niftiread(niiPath);
volume = double(squeeze(volume));

if ndims(volume) > 3
    volume = volume(:,:,:,1);
end

if ndims(volume) ~= 3
    error('El archivo NIfTI no es un volumen 3D válido.');
end

volume = volume - min(volume(:));

if max(volume(:)) > 0
    volume = volume / max(volume(:));
else
    error('Volumen sin intensidad útil.');
end

[N, M, Z] = size(volume);

if targetSlice < 1 || targetSlice > Z
    error('El corte %d no existe. El volumen tiene %d cortes.', targetSlice, Z);
end

%% Resolución espacial

if isfield(info, 'PixelDimensions') && length(info.PixelDimensions) >= 3
    dx = info.PixelDimensions(1);
    dy = info.PixelDimensions(2);
    dz = info.PixelDimensions(3);
else
    dx = NaN;
    dy = NaN;
    dz = NaN;
end

%% Extraer corte original

A = volume(:,:,targetSlice);
A = normalizar01(A);

%% ROI central

params.roiMinRow = 0.15;
params.roiMaxRow = 0.85;
params.roiMinCol = 0.15;
params.roiMaxCol = 0.85;

bodyMask = crearBodyMask(A);

ROI = crearROICentral(size(A), ...
    params.roiMinRow, params.roiMaxRow, ...
    params.roiMinCol, params.roiMaxCol);

ROI = ROI & bodyMask;

%% Banco de filtros

filtros = struct();

filtros.original = A;
filtros.mediana_3x3 = filtroMediana2D(A, 3);

kernelMedia = ones(3,3) / 9;
filtros.media_3x3 = filtroConvolucion2D(A, kernelMedia);

kernelG = kernelGaussiano(5, 1.0);
filtros.gaussiano_5x5 = filtroConvolucion2D(A, kernelG);

try
    filtros.wiener_5x5 = wiener2(A, [5 5]);
    filtros.wiener_5x5 = normalizar01(filtros.wiener_5x5);
catch
    filtros.wiener_5x5 = filtros.media_3x3;
end

baseSuave = filtros.gaussiano_5x5;
detalle = A - baseSuave;
filtros.unsharp_moderado = normalizar01(A + 0.5 * detalle);

%% Métricas de preprocesamiento

nombres = fieldnames(filtros);
n = length(nombres);

tabla = table();

for i = 1:n

    nombreFiltro = string(nombres{i});
    B = normalizar01(filtros.(nombres{i}));

    met = calcularMetricasPreprocesamiento(A, B, ROI);

    fila = table();

    fila.caseName = string(targetCase);
    fila.filename = string(niiPath);
    fila.slice = targetSlice;

    fila.dx_mm = dx;
    fila.dy_mm = dy;
    fila.dz_mm = dz;

    fila.filterName = nombreFiltro;

    fila.MSE_vs_original = met.mse;
    fila.PSNR_vs_original = met.psnr;
    fila.contrastROI = met.contrastROI;
    fila.noiseResidual = met.noiseResidual;
    fila.edgePreservation = met.edgePreservation;
    fila.sharpness = met.sharpness;
    fila.centralBrightCNR = met.centralBrightCNR;

    tabla = [tabla; fila];

end

contrastN = normalizarVector(tabla.contrastROI);
cnrN = normalizarVector(tabla.centralBrightCNR);
edgeN = normalizarVector(tabla.edgePreservation);
sharpN = normalizarVector(tabla.sharpness);
noiseN = normalizarVector(tabla.noiseResidual);
mseN = normalizarVector(tabla.MSE_vs_original);

tabla.preprocessingScore = ...
    0.25 * cnrN + ...
    0.20 * contrastN + ...
    0.20 * edgeN + ...
    0.15 * (1 - noiseN) + ...
    0.10 * sharpN + ...
    0.10 * (1 - mseN);

tabla.selectedForDetection = tabla.filterName == "mediana_3x3";
tabla = sortrows(tabla, 'preprocessingScore', 'descend');

%% Guardar imágenes

originalPath = fullfile(phase2Folder, ...
    'FaseII_original_t2_t25_00_corte93.png');

medianaPath = fullfile(phase2Folder, ...
    'FaseII_mediana3x3_t2_t25_00_corte93.png');

roiPath = fullfile(phase2Folder, ...
    'FaseII_ROI_central_t2_t25_00_corte93.png');

imwrite(uint8(255 * A), originalPath);
imwrite(uint8(255 * filtros.mediana_3x3), medianaPath);
imwrite(uint8(255 * ROI), roiPath);

%% Figura comparativa

fig1 = figure('Visible', 'on', ...
    'Position', [100 100 1700 850], ...
    'Color', 'w');

tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:n

    nombreFiltro = nombres{i};
    B = normalizar01(filtros.(nombreFiltro));

    nexttile;
    imshow(B, []);
    title(strrep(nombreFiltro, '_', ' '), 'Interpreter', 'none');

end

sgtitle(sprintf('FASE II | Preprocesamiento | %s | Corte %d', ...
    targetFile, targetSlice), ...
    'Interpreter', 'none', ...
    'FontSize', 15);

comparativaPath = fullfile(phase2Folder, ...
    'FaseII_comparativa_preprocesamiento_t2_t25_00_corte93.png');

saveas(fig1, comparativaPath);

%% Figura mediana

Bmed = filtros.mediana_3x3;
diffMed = abs(A - Bmed);

fig2 = figure('Visible', 'on', ...
    'Position', [100 100 1700 650], ...
    'Color', 'w');

subplot(1,4,1);
imshow(A, []);
title('Original corte 93', 'FontWeight', 'bold');

subplot(1,4,2);
imshow(Bmed, []);
title('Mediana 3x3', 'FontWeight', 'bold');

subplot(1,4,3);
imshow(diffMed, []);
title('Diferencia', 'FontWeight', 'bold');

subplot(1,4,4);
imshow(ROI, []);
title('ROI central', 'FontWeight', 'bold');

sgtitle('FASE II | Filtro de apoyo para detección: Mediana 3x3', ...
    'Interpreter', 'none', ...
    'FontSize', 15);

medianaResumenPath = fullfile(phase2Folder, ...
    'FaseII_resumen_mediana3x3_t2_t25_00_corte93.png');

saveas(fig2, medianaResumenPath);

%% CSV

csvPath = fullfile(phase2Folder, ...
    'FaseII_preprocesamiento_t2_t25_00_corte93.csv');

writetable(tabla, csvPath);

fprintf('\n====================================================\n');
fprintf('FASE 2 TERMINADA\n');
fprintf('Archivo: %s\n', targetFile);
fprintf('Corte: %d\n', targetSlice);
fprintf('Filtro seleccionado para apoyo: mediana_3x3\n');
fprintf('CSV:\n%s\n', csvPath);
fprintf('====================================================\n');

%% FUNCIONES LOCALES

function A = normalizar01(A)

    A = double(A);
    A = A - min(A(:));

    if max(A(:)) > 0
        A = A / max(A(:));
    end

end

function BW = crearBodyMask(A)

    Tbody = prctile(A(:), 10);
    BW = A > Tbody;

    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, round(0.005 * numel(A)));

    CC = bwconncomp(BW, 8);

    if CC.NumObjects > 0
        areas = cellfun(@numel, CC.PixelIdxList);
        [~, idx] = max(areas);

        BW2 = false(size(BW));
        BW2(CC.PixelIdxList{idx}) = true;
        BW = BW2;
    end

end

function ROI = crearROICentral(sz, rMin, rMax, cMin, cMax)

    N = sz(1);
    M = sz(2);

    r1 = max(1, round(rMin * N));
    r2 = min(N, round(rMax * N));

    c1 = max(1, round(cMin * M));
    c2 = min(M, round(cMax * M));

    ROIrect = false(N, M);
    ROIrect(r1:r2, c1:c2) = true;

    [X, Y] = meshgrid(1:M, 1:N);

    cx = (c1 + c2) / 2;
    cy = (r1 + r2) / 2;

    rx = (c2 - c1) / 2;
    ry = (r2 - r1) / 2;

    ROIellipse = ((X - cx).^2 / (rx^2 + eps)) + ...
                 ((Y - cy).^2 / (ry^2 + eps)) <= 1;

    ROI = ROIrect & ROIellipse;

end

function B = filtroMediana2D(A, tam)

    A = double(A);

    if mod(tam, 2) == 0
        tam = tam + 1;
    end

    [N, M] = size(A);
    p = floor(tam / 2);

    Ap = zeros(N + 2*p, M + 2*p);
    Ap(1+p:p+N, 1+p:p+M) = A;

    B = zeros(N, M);

    for i = 1:N
        for j = 1:M
            region = Ap(i:i+tam-1, j:j+tam-1);
            B(i,j) = median(region(:));
        end
    end

    B = normalizar01(B);

end

function B = filtroConvolucion2D(A, w)

    A = double(A);

    [N, M] = size(A);
    [nw, mw] = size(w);

    padN = floor(nw / 2);
    padM = floor(mw / 2);

    Ap = zeros(N + 2*padN, M + 2*padM);
    Ap(1+padN:padN+N, 1+padM:padM+M) = A;

    B = zeros(N, M);

    for i = 1:N
        for j = 1:M
            region = Ap(i:i+nw-1, j:j+mw-1);
            B(i,j) = sum(sum(region .* w));
        end
    end

    B = normalizar01(B);

end

function w = kernelGaussiano(tam, sigma)

    if mod(tam, 2) == 0
        tam = tam + 1;
    end

    centro = floor(tam / 2);
    w = zeros(tam, tam);

    for i = -centro:centro
        for j = -centro:centro
            w(i+centro+1, j+centro+1) = ...
                exp(-(i^2 + j^2) / (2 * sigma^2));
        end
    end

    w = w / sum(w(:));

end

function met = calcularMetricasPreprocesamiento(A, B, ROI)

    A = normalizar01(A);
    B = normalizar01(B);

    diff = A - B;

    met.mse = mean(diff(ROI).^2);
    met.psnr = 10 * log10(1 / (met.mse + eps));

    vals = B(ROI);
    met.contrastROI = std(vals);

    Bsuave = filtroConvolucion2D(B, ones(3,3)/9);
    residual = B - Bsuave;

    met.noiseResidual = std(residual(ROI));

    GA = magnitudGradiente(A);
    GB = magnitudGradiente(B);

    x = GA(ROI);
    y = GB(ROI);

    x = x - mean(x);
    y = y - mean(y);

    met.edgePreservation = sum(x .* y) / ...
        (sqrt(sum(x.^2) * sum(y.^2)) + eps);

    if isnan(met.edgePreservation)
        met.edgePreservation = 0;
    end

    met.sharpness = mean(GB(ROI));

    brightThr = prctile(vals, 88);
    brightVals = vals(vals >= brightThr);
    backVals = vals(vals < brightThr);

    if isempty(brightVals) || isempty(backVals)
        met.centralBrightCNR = 0;
    else
        met.centralBrightCNR = ...
            (mean(brightVals) - mean(backVals)) / (std(backVals) + eps);
    end

end

function G = magnitudGradiente(A)

    [Gx, Gy] = gradient(double(A));
    G = sqrt(Gx.^2 + Gy.^2);

end

function vNorm = normalizarVector(v)

    v = double(v);

    if max(v) - min(v) == 0
        vNorm = zeros(size(v));
    else
        vNorm = (v - min(v)) / (max(v) - min(v));
    end

end