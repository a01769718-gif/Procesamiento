clc;
clear;
close all;

%% ============================================================
% FASE 3 CORREGIDA
% Segmentación basada directamente en la lógica de FASE 1
%
% Caso: 5187149
% Corte: 77
%
% Cambio principal:
% - Se usa el corte original del recon.nii.
% - Se aplica filtro de mediana 2D manual antes de umbralizar.
% - Se usa la misma lógica de evaluarCorteVentricular del Código 1.
%
% Salida visual esperada:
% Original corte 77 | ROI central | Candidato ventricular | Overlay
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

phase3Folder = fullfile(rootFolder, ...
    'Resultados_Caso5187149_Corte77', ...
    'FaseIII_Segmentacion');

if ~exist(phase3Folder, 'dir')
    mkdir(phase3Folder);
end

%% Buscar recon.nii del caso 5187149

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

%% ============================================================
% Parámetros copiados de la lógica de Fase 1
%% ============================================================

params.minContentFrac = 0.04;
params.percentilBrillo = 88;
params.minAreaPx = 15;
params.maxAreaFracROI = 0.18;
params.scoreMin = 0.52;
params.topMaxPerCase = 8;
params.minSpacing = 2;

params.roiMinRow = 0.15;
params.roiMaxRow = 0.85;
params.roiMinCol = 0.15;
params.roiMaxCol = 0.85;

%% ============================================================
% Evaluar corte con la misma función de Fase 1
%% ============================================================

[met, BWvent, ROI] = evaluarCorteVentricular(A, params);

overlay = crearOverlayRojo(A, BWvent);

%% ============================================================
% Guardar imágenes individuales
%% ============================================================

originalPath = fullfile(phase3Folder, ...
    'FaseIII_original_caso5187149_corte77.png');

roiPath = fullfile(phase3Folder, ...
    'FaseIII_ROI_central_caso5187149_corte77.png');

maskPath = fullfile(phase3Folder, ...
    'FaseIII_mascara_caso5187149_corte77.png');

overlayPath = fullfile(phase3Folder, ...
    'FaseIII_overlay_caso5187149_corte77.png');

% Compatibilidad con Fase 4 anterior
compatImagePath = fullfile(phase3Folder, ...
    'FaseIII_imagen_filtrada_caso5187149_corte77.png');

imwrite(uint8(255 * A), originalPath);
imwrite(uint8(255 * A), compatImagePath);
imwrite(uint8(255 * ROI), roiPath);
imwrite(uint8(255 * BWvent), maskPath);
imwrite(overlay, overlayPath);

%% ============================================================
% Figura resumen como el Código 1
%% ============================================================

fig = figure('Visible', 'on', 'Position', [100 100 1500 600], 'Color', 'w');

subplot(1,4,1);
imshow(A, []);
title(['Original corte ', num2str(targetSlice)], 'FontWeight', 'bold');

subplot(1,4,2);
imshow(ROI, []);
title('ROI central', 'FontWeight', 'bold');

subplot(1,4,3);
imshow(BWvent, []);
title('Candidato ventricular', 'FontWeight', 'bold');

subplot(1,4,4);
imshow(overlay);
title('Overlay', 'FontWeight', 'bold');

sgtitle(['Caso ', targetCase, ...
    ' | Corte ', num2str(targetSlice), ...
    ' | VentricleScore = ', num2str(met.ventricleScore, '%.3f')], ...
    'Interpreter', 'none', ...
    'FontSize', 16);

summaryPath = fullfile(phase3Folder, ...
    'FaseIII_resumen_caso5187149_corte77.png');

saveas(fig, summaryPath);

%% ============================================================
% Guardar CSV
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

Tseg.contentFrac = met.contentFrac;
Tseg.candidateAreaPx = met.candidateAreaPx;
Tseg.areaFracROI = met.areaFracROI;
Tseg.numCandidates = met.numCandidates;
Tseg.CNR = met.CNR;
Tseg.centrality = met.centrality;
Tseg.symmetryScore = met.symmetryScore;
Tseg.shapeScore = met.shapeScore;
Tseg.ventricleScore = met.ventricleScore;

Tseg.originalPath = string(originalPath);
Tseg.roiPath = string(roiPath);
Tseg.maskPath = string(maskPath);
Tseg.overlayPath = string(overlayPath);
Tseg.summaryPath = string(summaryPath);

csvPath = fullfile(phase3Folder, ...
    'FaseIII_segmentacion_caso5187149_corte77.csv');

writetable(Tseg, csvPath);

fprintf('\n====================================================\n');
fprintf('FASE 3 TERMINADA\n');
fprintf('Caso: %s\n', targetCase);
fprintf('Corte: %d\n', targetSlice);
fprintf('VentricleScore: %.3f\n', met.ventricleScore);
fprintf('Área candidata: %d px\n', met.candidateAreaPx);
fprintf('Número de componentes: %d\n', met.numCandidates);
fprintf('Resumen visual:\n%s\n', summaryPath);
fprintf('Máscara:\n%s\n', maskPath);
fprintf('Overlay:\n%s\n', overlayPath);
fprintf('CSV:\n%s\n', csvPath);
fprintf('====================================================\n');

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

function [met, BWvalid, ROI] = evaluarCorteVentricular(A, params)

    A = normalizar01(A);

    [N, M] = size(A);

    %% Máscara de contenido anatómico

    Tbody = prctile(A(:), 10);
    bodyMask = A > Tbody;
    contentFrac = nnz(bodyMask) / numel(A);

    %% ROI central

    ROI = crearROICentral(size(A), ...
        params.roiMinRow, params.roiMaxRow, ...
        params.roiMinCol, params.roiMaxCol);

    ROI = ROI & bodyMask;

    if nnz(ROI) < 20
        BWvalid = false(size(A));
        met = metricasVacias(contentFrac);
        return;
    end

    %% ========================================================
    % FILTRO DE MEDIANA ANTES DE SELECCIONAR ÁREAS DE INTERÉS
    %
    % Esta es la parte clave tomada de Fase 1:
    % - A es el corte original.
    % - B es el corte suavizado con mediana.
    % - La detección del candidato se hace sobre B.
    % - El overlay final se dibuja sobre A.
    %% ========================================================

    B = filtroMediana2D(A, 3);

    valoresROI = B(ROI);

    T1 = prctile(valoresROI, params.percentilBrillo);
    T2 = mean(valoresROI) + 0.55 * std(valoresROI);

    threshold = max(T1, T2);

    BW = B > threshold;
    BW = BW & ROI;

    %% Limpieza morfológica

    BW = bwareaopen(BW, params.minAreaPx);

    try
        BW = imclose(BW, strel('disk', 2));
        BW = imopen(BW, strel('disk', 1));
    catch
        BW = imclose(BW, strel('square', 5));
        BW = imopen(BW, strel('square', 3));
    end

    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, params.minAreaPx);

    %% Componentes conectados

    CC = bwconncomp(BW, 8);

    if CC.NumObjects == 0
        BWvalid = false(size(A));
        met = metricasVacias(contentFrac);
        return;
    end

    stats = regionprops(CC, B, ...
        'Area', 'Centroid', 'MajorAxisLength', ...
        'MinorAxisLength', 'Eccentricity', 'MeanIntensity');

    maxAreaPx = params.maxAreaFracROI * nnz(ROI);

    centroImagen = [M/2, N/2];
    maxDist = sqrt((M/2)^2 + (N/2)^2);

    compScore = zeros(CC.NumObjects, 1);
    valido = false(CC.NumObjects, 1);

    outsideROI = B(ROI & ~BW);

    if isempty(outsideROI)
        outsideROI = B(ROI);
    end

    meanOutside = mean(outsideROI);
    stdOutside = std(outsideROI);

    for c = 1:CC.NumObjects

        area = stats(c).Area;

        if area < params.minAreaPx || area > maxAreaPx
            continue;
        end

        centroid = stats(c).Centroid;
        distancia = norm(centroid - centroImagen);
        centrality = 1 - distancia / (maxDist + eps);
        centrality = max(centrality, 0);

        meanInside = stats(c).MeanIntensity;
        cnrComp = (meanInside - meanOutside) / (stdOutside + eps);
        cnrScore = min(max(cnrComp / 3, 0), 1);

        ecc = stats(c).Eccentricity;
        shapeScore = min(max(ecc / 0.80, 0), 1);

        areaFracROI = area / (nnz(ROI) + eps);
        areaIdeal = 0.025;
        sigmaArea = 0.035;

        areaScore = exp(-((areaFracROI - areaIdeal)^2) / ...
            (2 * sigmaArea^2));

        compScore(c) = ...
            0.35 * cnrScore + ...
            0.25 * centrality + ...
            0.20 * areaScore + ...
            0.20 * shapeScore;

        valido(c) = true;

    end

    if ~any(valido)
        BWvalid = false(size(A));
        met = metricasVacias(contentFrac);
        return;
    end

    compScore(~valido) = -Inf;

    [~, orden] = sort(compScore, 'descend');

    maxComp = min(4, sum(valido));
    conservar = orden(1:maxComp);

    BWvalid = false(size(A));

    for k = 1:length(conservar)
        BWvalid(CC.PixelIdxList{conservar(k)}) = true;
    end

    BWvalid = BWvalid & ROI;

    %% Métricas globales

    candidateAreaPx = nnz(BWvalid);
    areaFracROI = candidateAreaPx / (nnz(ROI) + eps);

    CCvalid = bwconncomp(BWvalid, 8);
    numCandidates = CCvalid.NumObjects;

    inside = B(BWvalid);
    outside = B(ROI & ~BWvalid);

    if isempty(outside)
        outside = B(ROI);
    end

    if isempty(inside)
        CNR = 0;
    else
        CNR = (mean(inside) - mean(outside)) / (std(outside) + eps);
    end

    statsValid = regionprops(CCvalid, B, ...
        'Area', 'Centroid', 'Eccentricity', 'MeanIntensity');

    centrality = calcularCentralidad(statsValid, N, M);
    symmetryScore = calcularSimetria(statsValid, N, M);

    if isempty(statsValid)
        shapeScore = 0;
    else
        shapeScore = mean([statsValid.Eccentricity]);
    end

    %% Score final

    cnrScore = min(max(CNR / 3, 0), 1);

    areaIdeal = 0.045;
    sigmaArea = 0.050;

    areaScore = exp(-((areaFracROI - areaIdeal)^2) / ...
        (2 * sigmaArea^2));

    if numCandidates == 1 || numCandidates == 2
        componentScore = 1.0;
    elseif numCandidates <= 4
        componentScore = 0.75;
    elseif numCandidates <= 6
        componentScore = 0.40;
    else
        componentScore = 0.10;
    end

    ventricleScore = ...
        0.30 * cnrScore + ...
        0.20 * areaScore + ...
        0.20 * componentScore + ...
        0.15 * centrality + ...
        0.10 * symmetryScore + ...
        0.05 * shapeScore;

    ventricleScore = max(min(ventricleScore, 1), 0);

    %% Salida

    met.contentFrac = contentFrac;
    met.candidateAreaPx = candidateAreaPx;
    met.areaFracROI = areaFracROI;
    met.numCandidates = numCandidates;
    met.CNR = CNR;
    met.centrality = centrality;
    met.symmetryScore = symmetryScore;
    met.shapeScore = shapeScore;
    met.ventricleScore = ventricleScore;

end

function met = metricasVacias(contentFrac)

    met.contentFrac = contentFrac;
    met.candidateAreaPx = 0;
    met.areaFracROI = 0;
    met.numCandidates = 0;
    met.CNR = 0;
    met.centrality = 0;
    met.symmetryScore = 0;
    met.shapeScore = 0;
    met.ventricleScore = 0;

end

function centrality = calcularCentralidad(stats, N, M)

    if isempty(stats)
        centrality = 0;
        return;
    end

    centroImagen = [M/2, N/2];
    maxDist = sqrt((M/2)^2 + (N/2)^2);

    suma = 0;
    areaTotal = 0;

    for i = 1:length(stats)

        distancia = norm(stats(i).Centroid - centroImagen);
        c = 1 - distancia / (maxDist + eps);
        c = max(c, 0);

        suma = suma + c * stats(i).Area;
        areaTotal = areaTotal + stats(i).Area;

    end

    centrality = suma / (areaTotal + eps);

end

function symmetryScore = calcularSimetria(stats, N, M)

    if isempty(stats)
        symmetryScore = 0;
        return;
    end

    cx = arrayfun(@(s) s.Centroid(1), stats);
    cy = arrayfun(@(s) s.Centroid(2), stats);
    areas = arrayfun(@(s) s.Area, stats);

    leftIdx = find(cx < M/2);
    rightIdx = find(cx >= M/2);

    if ~isempty(leftIdx) && ~isempty(rightIdx)

        [~, il] = max(areas(leftIdx));
        [~, ir] = max(areas(rightIdx));

        L = leftIdx(il);
        R = rightIdx(ir);

        ySim = 1 - abs(cy(L) - cy(R)) / (0.30 * N + eps);
        ySim = max(min(ySim, 1), 0);

        areaSim = min(areas(L), areas(R)) / (max(areas(L), areas(R)) + eps);

        distL = abs(M/2 - cx(L));
        distR = abs(cx(R) - M/2);

        xBalance = 1 - abs(distL - distR) / (0.30 * M + eps);
        xBalance = max(min(xBalance, 1), 0);

        symmetryScore = 0.4 * ySim + 0.3 * areaSim + 0.3 * xBalance;

    elseif length(stats) == 1

        dCenter = abs(cx(1) - M/2) / (M/2 + eps);

        if dCenter < 0.15
            symmetryScore = 0.45;
        else
            symmetryScore = 0.20;
        end

    else

        symmetryScore = 0.25;

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

function overlay = crearOverlayRojo(A, BW)

    A = normalizar01(A);
    BW = logical(BW);

    borde = obtenerBordeBinario(BW);

    R = A;
    G = A;
    B = A;

    R(borde) = 1;
    G(borde) = 0;
    B(borde) = 0;

    interior = BW & ~borde;

    R(interior) = max(R(interior), 0.85);
    G(interior) = 0.5 * G(interior);
    B(interior) = 0.5 * B(interior);

    overlay = cat(3, R, G, B);
    overlay = uint8(255 * normalizar01(overlay));

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
        A = A / max(A(:));
    end

end