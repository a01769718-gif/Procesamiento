clc;
clear;
close all;

%% ============================================================
% RETO: MRI Fetal Image Processing for Ventriculomegaly Diagnosis
% CÓDIGO 1B: Selección de cortes con ventrículos notables
%
% Objetivo:
% - Revisar todos los recon.nii dentro de la carpeta principal.
% - Evaluar todos los cortes.
% - Seleccionar solo cortes con regiones brillantes centrales compatibles
%   con LCR/ventrículos laterales.
%
% IMPORTANTE:
% - Este código NO diagnostica ventriculomegalia.
% - Este código NO mide diámetro atrial.
% - Solo selecciona candidatos anatómicamente más prometedores.
%% ============================================================

%% 1. Seleccionar carpeta principal

rootFolder = uigetdir(pwd, ...
    'Selecciona la carpeta principal: placenta bad quality');

if rootFolder == 0
    error('No seleccionaste ninguna carpeta.');
end

%% 2. Crear carpeta de salida

outputFolder = fullfile(rootFolder, ...
    'Resultados_Cortes_Ventriculos_Notables');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% 3. Buscar archivos recon.nii

niiFiles = dir(fullfile(rootFolder, '**', 'recon.nii'));

if isempty(niiFiles)
    niiFiles = dir(fullfile(rootFolder, '**', '*.nii'));
end

if isempty(niiFiles)
    error('No se encontraron archivos .nii dentro de la carpeta seleccionada.');
end

fprintf('Se encontraron %d archivos NIfTI.\n', length(niiFiles));

%% 4. Parámetros de búsqueda ventricular

params.minContentFrac = 0.04;      % evita cortes casi vacíos
params.percentilBrillo = 88;       % busca regiones brillantes tipo LCR
params.minAreaPx = 15;             % elimina regiones muy pequeñas
params.maxAreaFracROI = 0.18;      % evita regiones demasiado grandes
params.scoreMin = 0.52;            % umbral para considerar ventrículo notable
params.topMaxPerCase = 8;          % máximo de cortes seleccionados por caso
params.minSpacing = 2;             % separación mínima entre cortes seleccionados

% ROI central aproximada
params.roiMinRow = 0.15;
params.roiMaxRow = 0.85;
params.roiMinCol = 0.15;
params.roiMaxCol = 0.85;

%% 5. Tablas de resultados

todosLosCortes = table();
cortesNotables = table();

%% 6. Procesar cada volumen

for fileIdx = 1:length(niiFiles)

    filename = fullfile(niiFiles(fileIdx).folder, niiFiles(fileIdx).name);
    [~, caseName] = fileparts(niiFiles(fileIdx).folder);

    fprintf('\n=============================================\n');
    fprintf('Procesando caso: %s\n', caseName);
    fprintf('Archivo: %s\n', filename);

    %% 6.1 Leer volumen

    try
        info = niftiinfo(filename);
        volume = niftiread(filename);
    catch
        warning('No se pudo leer el archivo: %s', filename);
        continue;
    end

    volume = double(squeeze(volume));

    if ndims(volume) > 3
        volume = volume(:,:,:,1);
    end

    if ndims(volume) ~= 3
        warning('El archivo no es un volumen 3D válido.');
        continue;
    end

    %% 6.2 Normalizar volumen completo

    volume = volume - min(volume(:));

    if max(volume(:)) > 0
        volume = volume / max(volume(:));
    else
        warning('Volumen sin intensidad útil.');
        continue;
    end

    [N, M, Z] = size(volume);

    if isfield(info, 'PixelDimensions') && length(info.PixelDimensions) >= 3
        dx = info.PixelDimensions(1);
        dy = info.PixelDimensions(2);
        dz = info.PixelDimensions(3);
    else
        dx = NaN;
        dy = NaN;
        dz = NaN;
    end

    %% 6.3 Evaluar todos los cortes

    metricasCaso = table();

    for z = 1:Z

        A = volume(:,:,z);
        A = normalizar01(A);

        [met, ~, ~] = evaluarCorteVentricular(A, params);

        if met.contentFrac < params.minContentFrac
            continue;
        end

        fila = table();

        fila.caseName = string(caseName);
        fila.filename = string(filename);
        fila.slice = z;

        fila.rows = N;
        fila.cols = M;
        fila.totalSlices = Z;

        fila.dx_mm = dx;
        fila.dy_mm = dy;
        fila.dz_mm = dz;

        fila.contentFrac = met.contentFrac;
        fila.candidateAreaPx = met.candidateAreaPx;
        fila.areaFracROI = met.areaFracROI;
        fila.numCandidates = met.numCandidates;
        fila.CNR = met.CNR;
        fila.centrality = met.centrality;
        fila.symmetryScore = met.symmetryScore;
        fila.shapeScore = met.shapeScore;
        fila.ventricleScore = met.ventricleScore;

        metricasCaso = [metricasCaso; fila];

    end

    if isempty(metricasCaso)
        warning('No se encontraron cortes útiles para el caso: %s', caseName);
        continue;
    end

    %% 6.4 Ordenar por score ventricular

    metricasCaso = sortrows(metricasCaso, 'ventricleScore', 'descend');

    %% 6.5 Seleccionar cortes con ventrículos notables

    seleccion = false(height(metricasCaso), 1);

    candidatos = find(metricasCaso.ventricleScore >= params.scoreMin);

    cortesYaElegidos = [];

    for k = 1:length(candidatos)

        idx = candidatos(k);
        corteActual = metricasCaso.slice(idx);

        if isempty(cortesYaElegidos)
            aceptar = true;
        else
            distanciaMin = min(abs(cortesYaElegidos - corteActual));
            aceptar = distanciaMin >= params.minSpacing;
        end

        if aceptar
            seleccion(idx) = true;
            cortesYaElegidos = [cortesYaElegidos; corteActual];
        end

        if sum(seleccion) >= params.topMaxPerCase
            break;
        end

    end

    metricasCaso.selectedVentricleCandidate = seleccion;

    seleccionadosCaso = metricasCaso(metricasCaso.selectedVentricleCandidate == true, :);

    fprintf('Cortes candidatos seleccionados: %d\n', height(seleccionadosCaso));

    %% 6.6 Crear carpeta del caso

    caseOutputFolder = fullfile(outputFolder, char(caseName));

    if ~exist(caseOutputFolder, 'dir')
        mkdir(caseOutputFolder);
    end

    %% 6.7 Guardar imágenes de cortes seleccionados

    for k = 1:height(seleccionadosCaso)

        z = seleccionadosCaso.slice(k);

        A = volume(:,:,z);
        A = normalizar01(A);

        [met, BWvent, ROI] = evaluarCorteVentricular(A, params);

        overlay = crearOverlayRojo(A, BWvent);

        fig = figure('Visible', 'off', 'Position', [100 100 1500 600]);

        subplot(1,4,1);
        imshow(A, []);
        title(['Original corte ', num2str(z)]);

        subplot(1,4,2);
        imshow(ROI, []);
        title('ROI central');

        subplot(1,4,3);
        imshow(BWvent, []);
        title('Candidato ventricular');

        subplot(1,4,4);
        imshow(overlay);
        title('Overlay');

        sgtitle(['Caso ', char(caseName), ...
            ' | Corte ', num2str(z), ...
            ' | VentricleScore = ', num2str(met.ventricleScore, '%.3f')], ...
            'Interpreter', 'none');

        saveName = fullfile(caseOutputFolder, ...
            ['Caso_', char(caseName), ...
            '_Corte_', num2str(z), ...
            '_VentriculoNotable_Score_', ...
            num2str(met.ventricleScore, '%.3f'), '.png']);

        saveas(fig, saveName);
        close(fig);

        imwrite(uint8(255 * A), fullfile(caseOutputFolder, ...
            ['Caso_', char(caseName), '_Corte_', num2str(z), '_original.png']));

        imwrite(uint8(255 * BWvent), fullfile(caseOutputFolder, ...
            ['Caso_', char(caseName), '_Corte_', num2str(z), '_mascara_candidata.png']));

        imwrite(overlay, fullfile(caseOutputFolder, ...
            ['Caso_', char(caseName), '_Corte_', num2str(z), '_overlay.png']));

    end

    %% 6.8 Guardar top 10 de revisión aunque no superen umbral

    nTopReview = min(10, height(metricasCaso));
    topReview = metricasCaso(1:nTopReview, :);

    writetable(topReview, fullfile(caseOutputFolder, ...
        ['top_revision_', char(caseName), '.csv']));

    %% 6.9 Agregar a tablas globales

    todosLosCortes = [todosLosCortes; metricasCaso];
    cortesNotables = [cortesNotables; seleccionadosCaso];

end

%% 7. Guardar CSV finales

csvTodos = fullfile(outputFolder, ...
    'todos_los_cortes_score_ventricular.csv');

csvNotables = fullfile(outputFolder, ...
    'cortes_ventriculos_notables.csv');

writetable(todosLosCortes, csvTodos);
writetable(cortesNotables, csvNotables);

fprintf('\n=============================================\n');
fprintf('Selección de cortes con ventrículos notables terminada.\n');
fprintf('Carpeta de salida:\n%s\n', outputFolder);
fprintf('Todos los cortes evaluados:\n%s\n', csvTodos);
fprintf('Cortes seleccionados:\n%s\n', csvNotables);
fprintf('=============================================\n');

disp('Cortes con ventrículos notables:');

if isempty(cortesNotables)
    disp('No se seleccionaron cortes automáticamente. Baja params.scoreMin o revisa top_revision por caso.');
else
    disp(cortesNotables(:, {'caseName', 'slice', ...
        'candidateAreaPx', 'numCandidates', ...
        'CNR', 'centrality', 'symmetryScore', ...
        'ventricleScore'}));
end

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

    %% Suavizado ligero para detección
    % No es la imagen final, solo ayuda a detectar regiones brillantes.

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

    %% Métricas globales de corte

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

    %% Score final del corte

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