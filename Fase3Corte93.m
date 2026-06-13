clc;
clear;
close all;

%% ============================================================
% FASE 3 CORREGIDA
% Segmentación automática anatómicamente guiada
%
% Archivo: t2-t25.00.nii
% Corte: 93
%
% Objetivo:
%   Detectar automáticamente los candidatos ventriculares laterales
%   usando filtro de mediana + ROI anatómica lateral.
%
% Importante:
%   NO se dibuja manualmente en esta fase.
%   La selección manual queda para Fase 4.
%
% Salida visual:
%   Original corte 93 | ROI ventricular | Candidato ventricular | Overlay
%% ============================================================

targetFile  = 't2-t25.00.nii';
targetCase  = 't2_t25_00';
targetSlice = 93;

%% Seleccionar carpeta

rootFolder = uigetdir(pwd, ...
    'Selecciona la carpeta MRI_t2w_nii');

if isequal(rootFolder, 0)
    error('No seleccionaste carpeta.');
end

phase3Folder = fullfile(rootFolder, ...
    'Resultados_t2_t25_00_Corte93', ...
    'FaseIII_Segmentacion');

if ~exist(phase3Folder, 'dir')
    mkdir(phase3Folder);
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
    volume = volume ./ max(volume(:));
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

%% ============================================================
% Parámetros ajustados para corte 93
%% ============================================================

params.percentilBrillo = 74;
params.stdFactor = 0.15;

params.minAreaPx = 30;
params.maxAreaFracROI = 0.45;

params.closeRadius = 3;
params.openRadius = 1;

params.maxComponentsPerROI = 1;

% ROI anatómica lateral.
% En este corte, los ventrículos reales están hacia la región lateral
% izquierda de la imagen, separados en una cavidad superior y una inferior.
params.supRowMin = 0.16;
params.supRowMax = 0.48;
params.supColMin = 0.06;
params.supColMax = 0.58;

params.infRowMin = 0.50;
params.infRowMax = 0.84;
params.infColMin = 0.06;
params.infColMax = 0.58;

% Se erosiona la máscara anatómica para evitar agarrar borde craneal.
params.bodyPercentile = 8;
params.bodyMinAreaFrac = 0.005;
params.bodyErodeRadius = 2;

%% ============================================================
% BodyMask
%% ============================================================

bodyMask = crearBodyMask(A, params);

try
    bodyInterior = imerode(bodyMask, strel('disk', params.bodyErodeRadius));
catch
    bodyInterior = imerode(bodyMask, strel('square', 2*params.bodyErodeRadius + 1));
end

if nnz(bodyInterior) < 50
    bodyInterior = bodyMask;
end

%% ============================================================
% ROI ventricular anatómica
%% ============================================================

ROI_sup = crearROIRectangular(size(A), ...
    params.supRowMin, params.supRowMax, ...
    params.supColMin, params.supColMax);

ROI_inf = crearROIRectangular(size(A), ...
    params.infRowMin, params.infRowMax, ...
    params.infColMin, params.infColMax);

ROI_sup = ROI_sup & bodyInterior;
ROI_inf = ROI_inf & bodyInterior;

ROI = ROI_sup | ROI_inf;

if nnz(ROI) < 30
    error('La ROI ventricular quedó muy pequeña. Ajusta los límites de ROI.');
end

%% ============================================================
% Filtro de mediana antes de seleccionar áreas de interés
%% ============================================================

B = filtroMediana2D(A, 3);

%% ============================================================
% Detectar ventrículo superior e inferior automáticamente
%% ============================================================

[BWsup, metSup] = detectarVentriculoEnROI(B, ROI_sup, params, 'superior');
[BWinf, metInf] = detectarVentriculoEnROI(B, ROI_inf, params, 'inferior');

BWvent = BWsup | BWinf;

BWvent = bwareaopen(BWvent, params.minAreaPx);

try
    BWvent = imclose(BWvent, strel('disk', 2));
catch
    BWvent = imclose(BWvent, strel('square', 5));
end

BWvent = imfill(BWvent, 'holes');
BWvent = BWvent & ROI;

%% ============================================================
% Métricas globales
%% ============================================================

CC = bwconncomp(BWvent, 8);
numCandidates = CC.NumObjects;

candidateAreaPx = nnz(BWvent);
areaFracROI = candidateAreaPx / (nnz(ROI) + eps);

inside = B(BWvent);
outside = B(ROI & ~BWvent);

if isempty(outside)
    outside = B(ROI);
end

if isempty(inside)
    CNR = 0;
else
    CNR = (mean(inside) - mean(outside)) / (std(outside) + eps);
end

statsValid = regionprops(CC, B, ...
    'Area', 'Centroid', 'Eccentricity', 'MeanIntensity');

centrality = calcularCentralidad(statsValid, N, M);
symmetryScore = calcularSimetriaPorDosRegiones(BWsup, BWinf, N, M);

if isempty(statsValid)
    shapeScore = 0;
else
    shapeScore = mean([statsValid.Eccentricity]);
end

cnrScore = min(max(CNR / 3, 0), 1);

areaIdeal = 0.18;
sigmaArea = 0.15;

areaScore = exp(-((areaFracROI - areaIdeal)^2) / ...
    (2 * sigmaArea^2));

if numCandidates == 2
    componentScore = 1.00;
elseif numCandidates == 1
    componentScore = 0.80;
elseif numCandidates <= 4
    componentScore = 0.55;
else
    componentScore = 0.20;
end

ventricleScore = ...
    0.30 * cnrScore + ...
    0.25 * areaScore + ...
    0.20 * componentScore + ...
    0.15 * centrality + ...
    0.10 * symmetryScore;

ventricleScore = max(min(ventricleScore, 1), 0);

%% ============================================================
% Overlay
%% ============================================================

overlay = crearOverlayDosColores(A, BWsup, BWinf);

%% ============================================================
% Guardar imágenes
%% ============================================================

originalPath = fullfile(phase3Folder, ...
    'FaseIII_original_t2_t25_00_corte93.png');

medianPath = fullfile(phase3Folder, ...
    'FaseIII_mediana3x3_t2_t25_00_corte93.png');

roiPath = fullfile(phase3Folder, ...
    'FaseIII_ROI_ventricular_t2_t25_00_corte93.png');

roiSupPath = fullfile(phase3Folder, ...
    'FaseIII_ROI_superior_t2_t25_00_corte93.png');

roiInfPath = fullfile(phase3Folder, ...
    'FaseIII_ROI_inferior_t2_t25_00_corte93.png');

maskSupPath = fullfile(phase3Folder, ...
    'FaseIII_mascara_ventriculo_superior_t2_t25_00_corte93.png');

maskInfPath = fullfile(phase3Folder, ...
    'FaseIII_mascara_ventriculo_inferior_t2_t25_00_corte93.png');

maskPath = fullfile(phase3Folder, ...
    'FaseIII_mascara_t2_t25_00_corte93.png');

overlayPath = fullfile(phase3Folder, ...
    'FaseIII_overlay_t2_t25_00_corte93.png');

imwrite(uint8(255 * A), originalPath);
imwrite(uint8(255 * B), medianPath);
imwrite(uint8(255 * ROI), roiPath);
imwrite(uint8(255 * ROI_sup), roiSupPath);
imwrite(uint8(255 * ROI_inf), roiInfPath);
imwrite(uint8(255 * BWsup), maskSupPath);
imwrite(uint8(255 * BWinf), maskInfPath);
imwrite(uint8(255 * BWvent), maskPath);
imwrite(overlay, overlayPath);

%% ============================================================
% Figura resumen estilo caso 77
%% ============================================================

summaryPath = fullfile(phase3Folder, ...
    'FaseIII_resumen_t2_t25_00_corte93.png');

fig = figure('Visible', 'on', ...
    'Position', [100 100 1700 650], ...
    'Color', 'w');

subplot(1,4,1);
imshow(A, []);
title('Original corte 93', 'FontWeight', 'bold');

subplot(1,4,2);
imshow(ROI, []);
title('ROI ventricular', 'FontWeight', 'bold');

subplot(1,4,3);
imshow(BWvent, []);
title('Candidato ventricular', 'FontWeight', 'bold');

subplot(1,4,4);
imshow(overlay);
title('Overlay', 'FontWeight', 'bold');

sgtitle(['Caso ', targetCase, ...
    ' | Corte ', num2str(targetSlice), ...
    ' | VentricleScore = ', num2str(ventricleScore, '%.3f')], ...
    'Interpreter', 'none', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

saveas(fig, summaryPath);

%% ============================================================
% CSV compatible con Fase 4
%% ============================================================

Tseg = table();

Tseg.caseName = string(targetCase);
Tseg.filename = string(niiPath);
Tseg.slice = targetSlice;

Tseg.rows = N;
Tseg.cols = M;
Tseg.totalSlices = Z;

Tseg.dx_mm = dx;
Tseg.dy_mm = dy;
Tseg.dz_mm = dz;

Tseg.segmentationMode = "automatico_roi_ventricular_lateral_mediana";

Tseg.areaSup_px = nnz(BWsup);
Tseg.areaInf_px = nnz(BWinf);
Tseg.candidateAreaPx = candidateAreaPx;
Tseg.areaFracROI = areaFracROI;

Tseg.numCandidates = numCandidates;
Tseg.CNR = CNR;
Tseg.centrality = centrality;
Tseg.symmetryScore = symmetryScore;
Tseg.shapeScore = shapeScore;
Tseg.ventricleScore = ventricleScore;

Tseg.scoreSup = metSup.score;
Tseg.scoreInf = metInf.score;

Tseg.originalPath = string(originalPath);
Tseg.medianPath = string(medianPath);
Tseg.roiPath = string(roiPath);
Tseg.roiSupPath = string(roiSupPath);
Tseg.roiInfPath = string(roiInfPath);
Tseg.maskSupPath = string(maskSupPath);
Tseg.maskInfPath = string(maskInfPath);
Tseg.maskPath = string(maskPath);
Tseg.overlayPath = string(overlayPath);
Tseg.summaryPath = string(summaryPath);

csvPath = fullfile(phase3Folder, ...
    'FaseIII_segmentacion_t2_t25_00_corte93.csv');

writetable(Tseg, csvPath);

fprintf('\n====================================================\n');
fprintf('FASE 3 CORREGIDA TERMINADA\n');
fprintf('Archivo: %s\n', targetFile);
fprintf('Corte: %d\n', targetSlice);
fprintf('Área superior: %d px\n', nnz(BWsup));
fprintf('Área inferior: %d px\n', nnz(BWinf));
fprintf('Área total candidata: %d px\n', candidateAreaPx);
fprintf('VentricleScore: %.3f\n', ventricleScore);
fprintf('CSV:\n%s\n', csvPath);
fprintf('Overlay:\n%s\n', overlayPath);
fprintf('Resumen:\n%s\n', summaryPath);
fprintf('====================================================\n');

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

function BW = crearBodyMask(A, params)

    Tbody = prctile(A(:), params.bodyPercentile);

    BW = A > Tbody;
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, round(params.bodyMinAreaFrac * numel(A)));

    CC = bwconncomp(BW, 8);

    if CC.NumObjects > 0

        areas = cellfun(@numel, CC.PixelIdxList);
        [~, idx] = max(areas);

        BW2 = false(size(BW));
        BW2(CC.PixelIdxList{idx}) = true;

        BW = BW2;

    end

end

function ROI = crearROIRectangular(sz, rMin, rMax, cMin, cMax)

    N = sz(1);
    M = sz(2);

    r1 = max(1, round(rMin * N));
    r2 = min(N, round(rMax * N));

    c1 = max(1, round(cMin * M));
    c2 = min(M, round(cMax * M));

    ROI = false(N, M);
    ROI(r1:r2, c1:c2) = true;

end

function [BWfinal, met] = detectarVentriculoEnROI(B, ROI, params, nombreROI)

    BWfinal = false(size(B));

    if nnz(ROI) < 20
        met.score = 0;
        met.area = 0;
        met.threshold = NaN;
        return;
    end

    valoresROI = B(ROI);

    T1 = prctile(valoresROI, params.percentilBrillo);
    T2 = mean(valoresROI) + params.stdFactor * std(valoresROI);

    threshold = max(T1, T2);

    BW = B > threshold;
    BW = BW & ROI;

    BW = bwareaopen(BW, params.minAreaPx);

    try
        BW = imclose(BW, strel('disk', params.closeRadius));
        BW = imopen(BW, strel('disk', params.openRadius));
    catch
        BW = imclose(BW, strel('square', 2*params.closeRadius + 1));
        BW = imopen(BW, strel('square', 2*params.openRadius + 1));
    end

    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, params.minAreaPx);

    CC = bwconncomp(BW, 8);

    if CC.NumObjects == 0
        met.score = 0;
        met.area = 0;
        met.threshold = threshold;
        return;
    end

    stats = regionprops(CC, B, ...
        'Area', 'Centroid', 'MeanIntensity', ...
        'Eccentricity', 'BoundingBox');

    maxAreaPx = params.maxAreaFracROI * nnz(ROI);

    outside = B(ROI & ~BW);

    if isempty(outside)
        outside = B(ROI);
    end

    meanOutside = mean(outside);
    stdOutside = std(outside);

    [rowsROI, colsROI] = find(ROI);
    centroROI = [mean(colsROI), mean(rowsROI)];

    maxDist = sqrt(size(B,1)^2 + size(B,2)^2);

    score = zeros(CC.NumObjects, 1);
    valido = false(CC.NumObjects, 1);

    for c = 1:CC.NumObjects

        areaC = stats(c).Area;

        if areaC < params.minAreaPx || areaC > maxAreaPx
            score(c) = -Inf;
            continue;
        end

        bbox = stats(c).BoundingBox;

        if bbox(3) < 4 || bbox(4) < 4
            score(c) = -Inf;
            continue;
        end

        centroid = stats(c).Centroid;

        distROI = norm(centroid - centroROI);
        centralityROI = 1 - distROI / (maxDist + eps);
        centralityROI = max(min(centralityROI, 1), 0);

        cnr = (stats(c).MeanIntensity - meanOutside) / (stdOutside + eps);
        cnrScore = min(max(cnr / 3, 0), 1);

        areaFrac = areaC / (nnz(ROI) + eps);

        areaScore = exp(-((areaFrac - 0.18)^2) / (2 * 0.12^2));

        eccScore = min(max(stats(c).Eccentricity / 0.95, 0), 1);

        score(c) = ...
            0.35 * cnrScore + ...
            0.30 * areaScore + ...
            0.20 * centralityROI + ...
            0.15 * eccScore;

        valido(c) = true;

    end

    if ~any(valido)
        met.score = 0;
        met.area = 0;
        met.threshold = threshold;
        return;
    end

    score(~valido) = -Inf;

    [~, orden] = sort(score, 'descend');

    keepN = min(params.maxComponentsPerROI, sum(valido));

    for k = 1:keepN
        idx = orden(k);
        BWfinal(CC.PixelIdxList{idx}) = true;
    end

    BWfinal = BWfinal & ROI;
    BWfinal = bwareaopen(BWfinal, params.minAreaPx);

    try
        BWfinal = imclose(BWfinal, strel('disk', 2));
    catch
        BWfinal = imclose(BWfinal, strel('square', 5));
    end

    BWfinal = imfill(BWfinal, 'holes');

    met.score = max(score(isfinite(score)));
    met.area = nnz(BWfinal);
    met.threshold = threshold;
    met.nombreROI = string(nombreROI);

end

function centrality = calcularCentralidad(stats, N, M)

    if isempty(stats)
        centrality = 0;
        return;
    end

    centro = [M/2, N/2];
    maxDist = sqrt((M/2)^2 + (N/2)^2);

    suma = 0;
    areaTotal = 0;

    for i = 1:length(stats)

        d = norm(stats(i).Centroid - centro);
        c = 1 - d / (maxDist + eps);
        c = max(c, 0);

        suma = suma + c * stats(i).Area;
        areaTotal = areaTotal + stats(i).Area;

    end

    centrality = suma / (areaTotal + eps);

end

function symmetryScore = calcularSimetriaPorDosRegiones(BWsup, BWinf, N, M)

    if nnz(BWsup) == 0 || nnz(BWinf) == 0
        symmetryScore = 0;
        return;
    end

    statsSup = regionprops(BWsup, 'Area', 'Centroid');
    statsInf = regionprops(BWinf, 'Area', 'Centroid');

    [~, idxS] = max([statsSup.Area]);
    [~, idxI] = max([statsInf.Area]);

    cS = statsSup(idxS).Centroid;
    cI = statsInf(idxI).Centroid;

    areaS = statsSup(idxS).Area;
    areaI = statsInf(idxI).Area;

    areaSim = min(areaS, areaI) / (max(areaS, areaI) + eps);

    xSim = 1 - abs(cS(1) - cI(1)) / (0.35 * M + eps);
    xSim = max(min(xSim, 1), 0);

    ySeparation = abs(cS(2) - cI(2)) / (0.35 * N + eps);
    ySeparation = max(min(ySeparation, 1), 0);

    symmetryScore = 0.40 * areaSim + 0.30 * xSim + 0.30 * ySeparation;
    symmetryScore = max(min(symmetryScore, 1), 0);

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

function overlay = crearOverlayDosColores(A, BWsup, BWinf)

    A = normalizar01(A);

    BWsup = logical(BWsup);
    BWinf = logical(BWinf);

    bordeSup = obtenerBordeBinario(BWsup);
    bordeInf = obtenerBordeBinario(BWinf);

    R = A;
    G = A;
    B = A;

    % Superior en amarillo
    R(bordeSup) = 1;
    G(bordeSup) = 1;
    B(bordeSup) = 0;

    intSup = BWsup & ~bordeSup;
    R(intSup) = max(R(intSup), 0.90);
    G(intSup) = max(G(intSup), 0.90);
    B(intSup) = 0.30 * B(intSup);

    % Inferior en cyan
    R(bordeInf) = 0;
    G(bordeInf) = 1;
    B(bordeInf) = 1;

    intInf = BWinf & ~bordeInf;
    R(intInf) = 0.30 * R(intInf);
    G(intInf) = max(G(intInf), 0.90);
    B(intInf) = max(B(intInf), 0.90);

    overlay = uint8(255 * cat(3, R, G, B));

end

function borde = obtenerBordeBinario(BW)

    BW = logical(BW);

    if nnz(BW) == 0
        borde = BW;
        return;
    end

    se = ones(3,3);
    vecinos = conv2(double(BW), se, 'same');

    erosionAprox = BW & (vecinos == 9);

    borde = BW & ~erosionAprox;

end

function A = normalizar01(A)

    A = double(A);
    A = A - min(A(:));

    if max(A(:)) > 0
        A = A ./ max(A(:));
    end

end