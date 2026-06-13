clc;
clear;
close all;

%% ============================================================
% FASE 2 CORREGIDA
% Preprocesamiento comparativo para:
% Caso: 5187149
% Corte: 77
%
% Cambio importante:
% - Esta fase YA NO entrega una imagen filtrada agresiva para Fase 3.
% - Solo documenta y compara filtros.
% - Se selecciona mediana 3x3 como filtro de apoyo para detección,
%   porque la Fase 3 usa el corte original + mediana interna.
%
% Salidas:
% - Figura comparativa de filtros
% - Figura enfocada en mediana 3x3
% - CSV con métricas
% - Imagen original y mediana guardadas para documentación
%% ============================================================

targetCase  = '5187149';
targetSlice = 77;

%% Seleccionar carpeta principal

rootFolder = uigetdir(pwd, ...
    'Selecciona la carpeta principal: placenta bad quality');

if isequal(rootFolder, 0)
    error('No seleccionaste carpeta.');
end

%% Carpeta de salida

phase2Folder = fullfile(rootFolder, ...
    'Resultados_Caso5187149_Corte77', ...
    'FaseII_Preprocesamiento');

if ~exist(phase2Folder, 'dir')
    mkdir(phase2Folder);
end

%% Buscar recon.nii del caso

niiFiles = dir(fullfile(rootFolder, '**', 'recon.nii'));

if isempty(niiFiles)
    niiFiles = dir(fullfile(rootFolder, '**', '*.nii'));
end

if isempty(niiFiles)
    error('No se encontraron archivos .nii dentro de la carpeta seleccionada.');
end

niiPath = "";

for k = 1:length(niiFiles)

    currentPath = fullfile(niiFiles(k).folder, niiFiles(k).name);

    if contains(string(currentPath), targetCase)
        niiPath = string(currentPath);
        break;
    end

end

if strlength(niiPath) == 0
    error('No se encontró recon.nii para el caso %s.', targetCase);
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

%% Crear ROI central igual a Fase 3

params.roiMinRow = 0.15;
params.roiMaxRow = 0.85;
params.roiMinCol = 0.15;
params.roiMaxCol = 0.85;

bodyMask = crearBodyMask(A);
ROI = crearROICentral(size(A), ...
    params.roiMinRow, params.roiMaxRow, ...
    params.roiMinCol, params.roiMaxCol);

ROI = ROI & bodyMask;

%% ============================================================
% Banco de filtros
%% ============================================================

filtros = struct();

filtros.original = A;

% Filtro de mediana manual, el mismo criterio que se usa en Fase 3.
filtros.mediana_3x3 = filtroMediana2D(A, 3);

% Filtro promedio 3x3.
kernelMedia = ones(3,3) / 9;
filtros.media_3x3 = filtroConvolucion2D(A, kernelMedia);

% Filtro gaussiano manual 5x5.
kernelG = kernelGaussiano(5, 1.0);
filtros.gaussiano_5x5 = filtroConvolucion2D(A, kernelG);

% Wiener si está disponible.
try
    filtros.wiener_5x5 = wiener2(A, [5 5]);
    filtros.wiener_5x5 = normalizar01(filtros.wiener_5x5);
catch
    filtros.wiener_5x5 = filtros.media_3x3;
end

% Realce unsharp moderado.
baseSuave = filtros.gaussiano_5x5;
detalle = A - baseSuave;
filtros.unsharp_moderado = normalizar01(A + 0.5 * detalle);

%% ============================================================
% Evaluar métricas de preprocesamiento
%% ============================================================

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

%% Score de utilidad para detección

% Normalización.
contrastN = normalizarVector(tabla.contrastROI);
cnrN = normalizarVector(tabla.centralBrightCNR);
edgeN = normalizarVector(tabla.edgePreservation);
sharpN = normalizarVector(tabla.sharpness);
noiseN = normalizarVector(tabla.noiseResidual);
mseN = normalizarVector(tabla.MSE_vs_original);

% En esta fase queremos:
% - conservar bordes,
% - mantener fidelidad al original,
% - reducir ruido,
% - mejorar contraste central,
% pero sin sobreprocesar.
tabla.preprocessingScore = ...
    0.25 * cnrN + ...
    0.20 * contrastN + ...
    0.20 * edgeN + ...
    0.15 * (1 - noiseN) + ...
    0.10 * sharpN + ...
    0.10 * (1 - mseN);

%% Forzar decisión metodológica

% Aunque se calculen todos los scores, para este pipeline se declara
% mediana_3x3 como filtro de apoyo, porque Fase 3 trabaja con:
% corte original + mediana 3x3 interna antes del umbral.
tabla.selectedForDetection = tabla.filterName == "mediana_3x3";

tabla = sortrows(tabla, 'preprocessingScore', 'descend');

%% Guardar imágenes individuales principales

originalPath = fullfile(phase2Folder, ...
    'FaseII_original_caso5187149_corte77.png');

medianaPath = fullfile(phase2Folder, ...
    'FaseII_mediana3x3_caso5187149_corte77.png');

roiPath = fullfile(phase2Folder, ...
    'FaseII_ROI_central_caso5187149_corte77.png');

imwrite(uint8(255 * A), originalPath);
imwrite(uint8(255 * filtros.mediana_3x3), medianaPath);
imwrite(uint8(255 * ROI), roiPath);

%% ============================================================
% Figura comparativa general
%% ============================================================

fig1 = figure('Visible', 'on', 'Position', [100 100 1700 850], 'Color', 'w');

tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:n

    nombreFiltro = nombres{i};
    B = normalizar01(filtros.(nombreFiltro));

    nexttile;
    imshow(B, []);
    title(strrep(nombreFiltro, '_', ' '), 'Interpreter', 'none');

end

sgtitle(sprintf('FASE II | Comparación de preprocesamiento | Caso %s | Corte %d', ...
    targetCase, targetSlice), ...
    'Interpreter', 'none', ...
    'FontSize', 15);

comparativaPath = fullfile(phase2Folder, ...
    'FaseII_comparativa_preprocesamiento_caso5187149_corte77.png');

saveas(fig1, comparativaPath);

%% ============================================================
% Figura enfocada en mediana 3x3
%% ============================================================

Bmed = filtros.mediana_3x3;
diffMed = abs(A - Bmed);

fig2 = figure('Visible', 'on', 'Position', [100 100 1700 650], 'Color', 'w');

subplot(1,4,1);
imshow(A, []);
title('Original corte 77', 'FontWeight', 'bold');

subplot(1,4,2);
imshow(Bmed, []);
title('Mediana 3x3', 'FontWeight', 'bold');

subplot(1,4,3);
imshow(diffMed, []);
title('Diferencia |Original - Mediana|', 'FontWeight', 'bold');

subplot(1,4,4);
imshow(ROI, []);
title('ROI central usada después', 'FontWeight', 'bold');

sgtitle(['FASE II | Filtro seleccionado para apoyo de detección: Mediana 3x3 | ', ...
    'Caso ', targetCase, ' | Corte ', num2str(targetSlice)], ...
    'Interpreter', 'none', ...
    'FontSize', 15);

medianaResumenPath = fullfile(phase2Folder, ...
    'FaseII_resumen_mediana3x3_caso5187149_corte77.png');

saveas(fig2, medianaResumenPath);

%% ============================================================
% Histogramas Original vs Mediana
%% ============================================================

fig3 = figure('Visible', 'on', 'Position', [100 100 1200 500], 'Color', 'w');

subplot(1,2,1);
histogram(A(ROI), 64);
title('Histograma ROI | Original');
xlim([0 1]);

subplot(1,2,2);
histogram(Bmed(ROI), 64);
title('Histograma ROI | Mediana 3x3');
xlim([0 1]);

sgtitle('FASE II | Comparación de histogramas dentro de ROI central', ...
    'Interpreter', 'none');

histPath = fullfile(phase2Folder, ...
    'FaseII_histograma_original_vs_mediana_caso5187149_corte77.png');

saveas(fig3, histPath);

%% ============================================================
% Guardar CSV
%% ============================================================

csvPath = fullfile(phase2Folder, ...
    'FaseII_preprocesamiento_caso5187149_corte77.csv');

writetable(tabla, csvPath);

%% Guardar resumen TXT

txtPath = fullfile(phase2Folder, ...
    'FaseII_decision_preprocesamiento_caso5187149_corte77.txt');

fid = fopen(txtPath, 'w');

fprintf(fid, 'FASE II - PREPROCESAMIENTO\n');
fprintf(fid, 'Caso: %s\n', targetCase);
fprintf(fid, 'Corte: %d\n\n', targetSlice);

fprintf(fid, 'Decision metodologica:\n');
fprintf(fid, 'Se selecciona mediana 3x3 como filtro de apoyo para deteccion.\n');
fprintf(fid, 'La imagen final de segmentacion no se reemplaza por la imagen filtrada.\n');
fprintf(fid, 'La Fase III usa el corte original y aplica mediana 3x3 internamente antes del umbral.\n\n');

fprintf(fid, 'Justificacion:\n');
fprintf(fid, '- Reduce ruido local sin deformar excesivamente la anatomia.\n');
fprintf(fid, '- Conserva mejor la visualizacion del corte axial.\n');
fprintf(fid, '- Evita el sobreprocesamiento observado con filtros mas agresivos.\n');
fprintf(fid, '- Mantiene coherencia con la logica de Fase I y Fase III.\n\n');

fprintf(fid, 'Archivos generados:\n');
fprintf(fid, 'Original: %s\n', originalPath);
fprintf(fid, 'Mediana: %s\n', medianaPath);
fprintf(fid, 'Comparativa: %s\n', comparativaPath);
fprintf(fid, 'Resumen mediana: %s\n', medianaResumenPath);
fprintf(fid, 'CSV: %s\n', csvPath);

fclose(fid);

%% Consola

fprintf('\n====================================================\n');
fprintf('FASE 2 CORREGIDA TERMINADA\n');
fprintf('Caso: %s\n', targetCase);
fprintf('Corte: %d\n', targetSlice);
fprintf('Filtro seleccionado para apoyo de detección: mediana_3x3\n');
fprintf('CSV:\n%s\n', csvPath);
fprintf('Comparativa:\n%s\n', comparativaPath);
fprintf('Resumen mediana:\n%s\n', medianaResumenPath);
fprintf('Resumen TXT:\n%s\n', txtPath);
fprintf('====================================================\n');

disp('Métricas de preprocesamiento:');
disp(tabla(:, {'filterName', 'MSE_vs_original', 'contrastROI', ...
    'noiseResidual', 'edgePreservation', 'centralBrightCNR', ...
    'preprocessingScore', 'selectedForDetection'}));

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

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