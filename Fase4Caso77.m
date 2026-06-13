clc;
clear;
close all;

%% ============================================================
% FASE 4 CORREGIDA
% Cuantificación manual/semimanual de DOS ventrículos
%
% Caso: 5187149
% Corte: 77
%
% Uso:
% 1. Selecciona la carpeta "placenta bad quality".
% 2. Dibuja el ventrículo 1.
% 3. Dibuja el ventrículo 2.
% 4. Dibuja el diámetro atrial del ventrículo 1.
% 5. Dibuja el diámetro atrial del ventrículo 2.
%
% Calcula:
% - Diámetro atrial de cada ventrículo en px y mm.
% - Diámetro atrial máximo.
% - Área de cada ventrículo.
% - Área total ventricular.
% - Volumen estimado del corte.
% - Eje mayor/eje menor.
% - Excentricidad.
% - Relación eje mayor/eje menor.
% - Compacidad.
% - Clasificación clínica preliminar.
%% ============================================================

targetCase  = '5187149';
targetSlice = 77;

%% Seleccionar carpeta principal

rootFolder = uigetdir(pwd, ...
    'Selecciona la carpeta principal: placenta bad quality');

if isequal(rootFolder, 0)
    error('No seleccionaste carpeta.');
end

%% Carpetas

phase3Folder = fullfile(rootFolder, ...
    'Resultados_Caso5187149_Corte77', ...
    'FaseIII_Segmentacion');

phase4Folder = fullfile(rootFolder, ...
    'Resultados_Caso5187149_Corte77', ...
    'FaseIV_Cuantificacion_Manual_Dos_Ventriculos');

if ~exist(phase4Folder, 'dir')
    mkdir(phase4Folder);
end

%% Archivos de entrada

csvFase3 = fullfile(phase3Folder, ...
    'FaseIII_segmentacion_caso5187149_corte77.csv');

originalPath = fullfile(phase3Folder, ...
    'FaseIII_original_caso5187149_corte77.png');

if ~exist(csvFase3, 'file')
    error('No se encontró el CSV de Fase 3:\n%s', csvFase3);
end

if ~exist(originalPath, 'file')
    error('No se encontró la imagen original de Fase 3:\n%s', originalPath);
end

%% Leer datos de Fase 3

T3 = readtable(csvFase3, 'TextType', 'string');

dx = T3.dx_mm(1);
dy = T3.dy_mm(1);
dz = T3.dz_mm(1);

if isnan(dx) || isnan(dy) || dx <= 0 || dy <= 0
    warning('No hay resolución espacial válida. Se usará 1 mm/px.');
    dx = 1;
    dy = 1;
end

if isnan(dz) || dz <= 0
    warning('No hay espesor de corte válido. Se usará dz = 1 mm para volumen estimado.');
    dz = 1;
end

pixelSizeMean = mean([dx, dy]);

%% Leer imagen original

I = imread(originalPath);

if ndims(I) == 3
    I = rgb2gray(I);
end

I = im2double(I);
I = normalizar01(I);

%% ============================================================
% PASO 1: Dibujar los dos ventrículos
%% ============================================================

ventNames = ["Ventriculo_1"; "Ventriculo_2"];
BWvent = cell(2,1);

for v = 1:2

    figMask = figure('Name', ['Dibuja ', char(ventNames(v))], ...
        'Position', [100 100 950 720], ...
        'Color', 'w');

    imshow(I, []);
    hold on;

    if v == 2 && ~isempty(BWvent{1})
        visboundaries(BWvent{1}, 'Color', 'y', 'LineWidth', 1.5);
    end

    title({['Dibuja ', char(ventNames(v))], ...
           'Doble clic para terminar el trazo'}, ...
           'FontSize', 13, 'FontWeight', 'bold');

    h = drawfreehand('Color', 'yellow', 'LineWidth', 2);
    BW = createMask(h);

    close(figMask);

    BW = logical(BW);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, 3);

    if nnz(BW) == 0
        error('La máscara de %s quedó vacía. Vuelve a correr y dibújala de nuevo.', ventNames(v));
    end

    BWvent{v} = BW;

end

BW1 = BWvent{1};
BW2 = BWvent{2};
BWtotal = BW1 | BW2;

%% ============================================================
% PASO 2: Dibujar diámetro atrial para cada ventrículo
%% ============================================================

lineas = zeros(2,4);
diametro_px = zeros(2,1);
diametro_mm = zeros(2,1);

for v = 1:2

    figLine = figure('Name', ['Dibuja diámetro atrial ', char(ventNames(v))], ...
        'Position', [100 100 950 720], ...
        'Color', 'w');

    imshow(I, []);
    hold on;

    visboundaries(BWvent{v}, 'Color', 'y', 'LineWidth', 1.5);

    title({['Dibuja el diámetro atrial de ', char(ventNames(v))], ...
           'Traza la línea de borde a borde del ventrículo'}, ...
           'FontSize', 13, 'FontWeight', 'bold');

    hLine = drawline('Color', 'cyan', 'LineWidth', 2);

    pos = hLine.Position;

    x1 = pos(1,1);
    y1 = pos(1,2);
    x2 = pos(2,1);
    y2 = pos(2,2);

    lineas(v,:) = [x1 y1 x2 y2];

    diametro_px(v) = sqrt((x2 - x1)^2 + (y2 - y1)^2);

    diametro_mm(v) = sqrt(((x2 - x1) * dx)^2 + ((y2 - y1) * dy)^2);

    close(figLine);

end

%% ============================================================
% Métricas por ventrículo
%% ============================================================

tablaVentriculos = table();

for v = 1:2

    BW = BWvent{v};

    CC = bwconncomp(BW, 8);

    stats = regionprops(CC, ...
        'Area', ...
        'Perimeter', ...
        'MajorAxisLength', ...
        'MinorAxisLength', ...
        'Eccentricity', ...
        'Orientation', ...
        'Centroid', ...
        'BoundingBox');

    if CC.NumObjects == 0
        error('No se detectó componente en %s.', ventNames(v));
    end

    areas = [stats.Area];
    [~, idxPrincipal] = max(areas);

    S = stats(idxPrincipal);

    area_px = nnz(BW);
    area_mm2 = area_px * dx * dy;

    perimeter_px = 0;

    for k = 1:length(stats)
        perimeter_px = perimeter_px + stats(k).Perimeter;
    end

    majorAxis_px = S.MajorAxisLength;
    minorAxis_px = S.MinorAxisLength;

    majorAxis_mm = majorAxis_px * pixelSizeMean;
    minorAxis_mm = minorAxis_px * pixelSizeMean;

    eccentricity = S.Eccentricity;

    relacionMayorMenor = majorAxis_mm / (minorAxis_mm + eps);

    compactness = (4 * pi * area_px) / (perimeter_px^2 + eps);

    diametroEquivalente_mm = 2 * sqrt(area_mm2 / pi);

    volumenEstimado_mm3 = area_mm2 * dz;

    fila = table();

    fila.caseName = string(targetCase);
    fila.slice = targetSlice;
    fila.ventriculo = ventNames(v);

    fila.dx_mm = dx;
    fila.dy_mm = dy;
    fila.dz_mm = dz;

    fila.diametroAtrial_px = diametro_px(v);
    fila.diametroAtrial_mm = diametro_mm(v);

    fila.area_px = area_px;
    fila.area_mm2 = area_mm2;
    fila.volumenEstimadoCorte_mm3 = volumenEstimado_mm3;

    fila.perimeter_px = perimeter_px;

    fila.majorAxis_px = majorAxis_px;
    fila.minorAxis_px = minorAxis_px;

    fila.majorAxis_mm = majorAxis_mm;
    fila.minorAxis_mm = minorAxis_mm;

    fila.relacionEjeMayorMenor = relacionMayorMenor;
    fila.eccentricity = eccentricity;
    fila.compactness = compactness;
    fila.diametroEquivalente_mm = diametroEquivalente_mm;

    fila.centroidX = S.Centroid(1);
    fila.centroidY = S.Centroid(2);
    fila.orientation_deg = S.Orientation;

    tablaVentriculos = [tablaVentriculos; fila];

end

%% ============================================================
% Métricas globales de ambos ventrículos
%% ============================================================

areaTotal_px = nnz(BWtotal);
areaTotal_mm2 = areaTotal_px * dx * dy;
volumenTotalEstimado_mm3 = areaTotal_mm2 * dz;

diametroAtrialMax_mm = max(diametro_mm);
diametroAtrialPromedio_mm = mean(diametro_mm);

[~, idxMaxDiam] = max(diametro_mm);
ventriculoCritico = ventNames(idxMaxDiam);

diagnostico = clasificarVentriculomegalia(diametroAtrialMax_mm);

if diametroAtrialMax_mm < 10
    riesgo = "Sin sospecha preliminar";
else
    riesgo = "Sospecha por criterio >= 10 mm";
end

eccentricityMean = mean(tablaVentriculos.eccentricity, 'omitnan');
relacionMean = mean(tablaVentriculos.relacionEjeMayorMenor, 'omitnan');
compactnessMean = mean(tablaVentriculos.compactness, 'omitnan');

%% Interpretación morfológica

if eccentricityMean < 0.50
    interpretacionEcc = "Morfología global relativamente redondeada.";
elseif eccentricityMean < 0.80
    interpretacionEcc = "Morfología global moderadamente elíptica.";
else
    interpretacionEcc = "Morfología global alargada.";
end

if relacionMean < 1.5
    interpretacionRelacion = "Relación eje mayor/eje menor baja; cavidades poco alargadas.";
elseif relacionMean < 2.5
    interpretacionRelacion = "Relación eje mayor/eje menor moderada; cavidades elípticas.";
else
    interpretacionRelacion = "Relación eje mayor/eje menor alta; cavidades muy alargadas.";
end

if compactnessMean > 0.70
    interpretacionCompacidad = "Alta compacidad promedio.";
elseif compactnessMean > 0.35
    interpretacionCompacidad = "Compacidad intermedia; forma elíptica o algo irregular.";
else
    interpretacionCompacidad = "Baja compacidad; posible fragmentación o irregularidad.";
end

%% ============================================================
% Guardar máscaras
%% ============================================================

maskV1Path = fullfile(phase4Folder, ...
    'FaseIV_mascara_manual_ventriculo1_caso5187149_corte77.png');

maskV2Path = fullfile(phase4Folder, ...
    'FaseIV_mascara_manual_ventriculo2_caso5187149_corte77.png');

maskTotalPath = fullfile(phase4Folder, ...
    'FaseIV_mascara_manual_dos_ventriculos_caso5187149_corte77.png');

imwrite(uint8(255 * BW1), maskV1Path);
imwrite(uint8(255 * BW2), maskV2Path);
imwrite(uint8(255 * BWtotal), maskTotalPath);

%% ============================================================
% Visual principal
%% ============================================================

visualPath = fullfile(phase4Folder, ...
    'FaseIV_visual_cuantificacion_dos_ventriculos_caso5187149_corte77.png');

fig = figure('Visible', 'on', ...
    'Position', [100 100 1900 900], ...
    'Color', 'w');

subplot(2,3,1);
imshow(I, []);
title('Original corte 77', 'FontWeight', 'bold');

subplot(2,3,2);
imshow(BWtotal, []);
title('Máscara manual: ambos ventrículos', 'FontWeight', 'bold');

subplot(2,3,3);
imshow(I, []);
hold on;
visboundaries(BW1, 'Color', 'y', 'LineWidth', 1.5);
visboundaries(BW2, 'Color', 'c', 'LineWidth', 1.5);
title('Overlay ambos ventrículos', 'FontWeight', 'bold');
hold off;

subplot(2,3,4);
imshow(I, []);
hold on;
visboundaries(BW1, 'Color', 'y', 'LineWidth', 1.5);
plot([lineas(1,1) lineas(1,3)], [lineas(1,2) lineas(1,4)], ...
    'c-', 'LineWidth', 3);
text(mean([lineas(1,1), lineas(1,3)]) + 5, ...
     mean([lineas(1,2), lineas(1,4)]), ...
     [num2str(diametro_mm(1), '%.2f'), ' mm'], ...
     'Color', 'cyan', ...
     'FontSize', 12, ...
     'FontWeight', 'bold', ...
     'BackgroundColor', 'black');
title('Diámetro atrial ventrículo 1', 'FontWeight', 'bold');
hold off;

subplot(2,3,5);
imshow(I, []);
hold on;
visboundaries(BW2, 'Color', 'c', 'LineWidth', 1.5);
plot([lineas(2,1) lineas(2,3)], [lineas(2,2) lineas(2,4)], ...
    'y-', 'LineWidth', 3);
text(mean([lineas(2,1), lineas(2,3)]) + 5, ...
     mean([lineas(2,2), lineas(2,4)]), ...
     [num2str(diametro_mm(2), '%.2f'), ' mm'], ...
     'Color', 'yellow', ...
     'FontSize', 12, ...
     'FontWeight', 'bold', ...
     'BackgroundColor', 'black');
title('Diámetro atrial ventrículo 2', 'FontWeight', 'bold');
hold off;

subplot(2,3,6);
axis off;

texto = {
    ['Caso: ', targetCase]
    ['Corte: ', num2str(targetSlice)]
    ''
    ['Pixel: ', num2str(dx, '%.4f'), ' × ', num2str(dy, '%.4f'), ' mm']
    ['Área por pixel: ', num2str(dx*dy, '%.4f'), ' mm²']
    ''
    ['Diámetro V1: ', num2str(diametro_mm(1), '%.2f'), ' mm']
    ['Diámetro V2: ', num2str(diametro_mm(2), '%.2f'), ' mm']
    ['Diámetro máximo: ', num2str(diametroAtrialMax_mm, '%.2f'), ' mm']
    ['Ventrículo crítico: ', char(ventriculoCritico)]
    ''
    ['Diagnóstico preliminar: ', char(diagnostico)]
    ['Riesgo: ', char(riesgo)]
    ''
    ['Área total: ', num2str(areaTotal_px), ' px']
    ['Área total: ', num2str(areaTotal_mm2, '%.2f'), ' mm²']
    ['Volumen estimado corte: ', num2str(volumenTotalEstimado_mm3, '%.2f'), ' mm³']
    ''
    ['Excentricidad media: ', num2str(eccentricityMean, '%.3f')]
    ['Relación eje mayor/eje menor media: ', num2str(relacionMean, '%.3f')]
    ['Compacidad media: ', num2str(compactnessMean, '%.3f')]
};

text(0.02, 0.98, texto, ...
    'VerticalAlignment', 'top', ...
    'FontSize', 11.5, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none');

sgtitle(['FASE IV | Cuantificación manual de dos ventrículos | Caso ', targetCase, ...
    ' | Corte ', num2str(targetSlice)], ...
    'Interpreter', 'none', ...
    'FontSize', 16);

saveas(fig, visualPath);

%% ============================================================
% Visual de conversión px a mm
%% ============================================================

conversionPath = fullfile(phase4Folder, ...
    'FaseIV_conversion_px_mm_dos_ventriculos_caso5187149_corte77.png');

fig2 = figure('Visible', 'on', ...
    'Position', [100 100 1600 650], ...
    'Color', 'w');

subplot(1,3,1);
imshow(I, []);
title('Original');

subplot(1,3,2);
imshow(BWtotal, []);
title('Máscara ambos ventrículos');

subplot(1,3,3);
imshow(I, []);
hold on;
visboundaries(BW1, 'Color', 'y', 'LineWidth', 1.5);
visboundaries(BW2, 'Color', 'c', 'LineWidth', 1.5);

plot([lineas(1,1) lineas(1,3)], [lineas(1,2) lineas(1,4)], ...
    'c-', 'LineWidth', 3);
plot([lineas(2,1) lineas(2,3)], [lineas(2,2) lineas(2,4)], ...
    'y-', 'LineWidth', 3);

barra_mm = 5;
barra_px = barra_mm / pixelSizeMean;

x0 = round(0.08 * size(I,2));
y0 = round(0.90 * size(I,1));

plot([x0 x0 + barra_px], [y0 y0], 'w-', 'LineWidth', 4);
text(x0, y0 - 8, [num2str(barra_mm), ' mm'], ...
    'Color', 'white', ...
    'FontWeight', 'bold', ...
    'BackgroundColor', 'black');

text(10, 20, {
    ['dx = ', num2str(dx, '%.4f'), ' mm/px']
    ['dy = ', num2str(dy, '%.4f'), ' mm/px']
    ['Área/pixel = ', num2str(dx*dy, '%.4f'), ' mm²']
    ['Diámetro V1 = ', num2str(diametro_mm(1), '%.2f'), ' mm']
    ['Diámetro V2 = ', num2str(diametro_mm(2), '%.2f'), ' mm']
    ['Diámetro máximo = ', num2str(diametroAtrialMax_mm, '%.2f'), ' mm']
    ['Área total = ', num2str(areaTotal_mm2, '%.2f'), ' mm²']
    }, ...
    'Color', 'yellow', ...
    'FontWeight', 'bold', ...
    'BackgroundColor', 'black', ...
    'Interpreter', 'none');

title('Conversión px → mm');
hold off;

sgtitle('FASE IV | Conversión espacial de ambos ventrículos', ...
    'Interpreter', 'none', ...
    'FontSize', 15);

saveas(fig2, conversionPath);

%% ============================================================
% Guardar CSV por ventrículo
%% ============================================================

csvVentriculos = fullfile(phase4Folder, ...
    'FaseIV_metricas_por_ventriculo_caso5187149_corte77.csv');

writetable(tablaVentriculos, csvVentriculos);

%% Guardar CSV resumen
%% ============================================================

Tresumen = table();

Tresumen.caseName = string(targetCase);
Tresumen.filename = T3.filename(1);
Tresumen.slice = targetSlice;

Tresumen.dx_mm = dx;
Tresumen.dy_mm = dy;
Tresumen.dz_mm = dz;

Tresumen.diametroVentriculo1_mm = diametro_mm(1);
Tresumen.diametroVentriculo2_mm = diametro_mm(2);
Tresumen.diametroAtrialMax_mm = diametroAtrialMax_mm;
Tresumen.diametroAtrialPromedio_mm = diametroAtrialPromedio_mm;
Tresumen.ventriculoCritico = ventriculoCritico;

Tresumen.diagnosticoPreliminar = diagnostico;
Tresumen.riesgo = riesgo;

Tresumen.areaTotal_px = areaTotal_px;
Tresumen.areaTotal_mm2 = areaTotal_mm2;
Tresumen.volumenTotalEstimadoCorte_mm3 = volumenTotalEstimado_mm3;

Tresumen.eccentricityMean = eccentricityMean;
Tresumen.relacionEjeMayorMenorMean = relacionMean;
Tresumen.compactnessMean = compactnessMean;

Tresumen.interpretacionEccentricity = interpretacionEcc;
Tresumen.interpretacionRelacionEjes = interpretacionRelacion;
Tresumen.interpretacionCompacidad = interpretacionCompacidad;

Tresumen.maskV1Path = string(maskV1Path);
Tresumen.maskV2Path = string(maskV2Path);
Tresumen.maskTotalPath = string(maskTotalPath);
Tresumen.visualCuantificacionPath = string(visualPath);
Tresumen.visualConversionPath = string(conversionPath);

csvResumen = fullfile(phase4Folder, ...
    'FaseIV_resumen_dos_ventriculos_caso5187149_corte77.csv');

writetable(Tresumen, csvResumen);

%% ============================================================
% Guardar TXT
%% ============================================================

txtPath = fullfile(phase4Folder, ...
    'FaseIV_resumen_dos_ventriculos_caso5187149_corte77.txt');

fid = fopen(txtPath, 'w');

fprintf(fid, 'FASE IV - CUANTIFICACION MANUAL DE DOS VENTRICULOS\n');
fprintf(fid, 'Caso: %s\n', targetCase);
fprintf(fid, 'Corte: %d\n\n', targetSlice);

fprintf(fid, 'Resolucion espacial:\n');
fprintf(fid, 'dx = %.6f mm/px\n', dx);
fprintf(fid, 'dy = %.6f mm/px\n', dy);
fprintf(fid, 'dz = %.6f mm\n', dz);
fprintf(fid, 'Area por pixel = %.6f mm2\n\n', dx*dy);

fprintf(fid, 'Diametros atriales:\n');
fprintf(fid, 'Ventriculo 1 = %.4f mm\n', diametro_mm(1));
fprintf(fid, 'Ventriculo 2 = %.4f mm\n', diametro_mm(2));
fprintf(fid, 'Diametro maximo = %.4f mm\n', diametroAtrialMax_mm);
fprintf(fid, 'Ventriculo critico = %s\n', ventriculoCritico);
fprintf(fid, 'Diagnostico preliminar = %s\n', diagnostico);
fprintf(fid, 'Riesgo = %s\n\n', riesgo);

fprintf(fid, 'Area y volumen:\n');
fprintf(fid, 'Area total = %d px\n', areaTotal_px);
fprintf(fid, 'Area total = %.4f mm2\n', areaTotal_mm2);
fprintf(fid, 'Volumen estimado del corte = %.4f mm3\n', volumenTotalEstimado_mm3);
fprintf(fid, 'Nota: al tener un solo corte, este volumen es parcial, no volumen ventricular total.\n\n');

fprintf(fid, 'Descriptores morfologicos promedio:\n');
fprintf(fid, 'Excentricidad media = %.4f\n', eccentricityMean);
fprintf(fid, 'Relacion eje mayor/eje menor media = %.4f\n', relacionMean);
fprintf(fid, 'Compacidad media = %.4f\n\n', compactnessMean);

fprintf(fid, 'Interpretacion:\n');
fprintf(fid, '%s\n', interpretacionEcc);
fprintf(fid, '%s\n', interpretacionRelacion);
fprintf(fid, '%s\n\n', interpretacionCompacidad);

fprintf(fid, 'Relacion con severidad:\n');
fprintf(fid, ['Al aumentar la severidad de la ventriculomegalia, ', ...
    'se espera incremento del diametro atrial y del area ventricular. ', ...
    'Los descriptores morfologicos pueden cambiar si la dilatacion produce ', ...
    'cavidades mas elongadas o irregulares, aumentando la relacion eje mayor/eje menor ', ...
    'y modificando excentricidad y compacidad.\n\n']);

fprintf(fid, 'Advertencia:\n');
fprintf(fid, 'Resultado manual/semiautomatico preliminar. No sustituye diagnostico medico.\n');

fclose(fid);

%% Consola

fprintf('\n====================================================\n');
fprintf('FASE 4 DOS VENTRICULOS TERMINADA\n');
fprintf('Caso: %s\n', targetCase);
fprintf('Corte: %d\n', targetSlice);
fprintf('Diametro V1: %.2f mm\n', diametro_mm(1));
fprintf('Diametro V2: %.2f mm\n', diametro_mm(2));
fprintf('Diametro maximo: %.2f mm\n', diametroAtrialMax_mm);
fprintf('Diagnostico preliminar: %s\n', diagnostico);
fprintf('Area total: %.2f mm2\n', areaTotal_mm2);
fprintf('Volumen estimado del corte: %.2f mm3\n', volumenTotalEstimado_mm3);
fprintf('CSV por ventriculo:\n%s\n', csvVentriculos);
fprintf('CSV resumen:\n%s\n', csvResumen);
fprintf('Visual cuantificacion:\n%s\n', visualPath);
fprintf('Visual conversion:\n%s\n', conversionPath);
fprintf('Resumen TXT:\n%s\n', txtPath);
fprintf('====================================================\n');

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

function diagnostico = clasificarVentriculomegalia(d)

    if isnan(d) || d <= 0
        diagnostico = "No cuantificable";
    elseif d < 10
        diagnostico = "Normal preliminar";
    elseif d <= 12
        diagnostico = "Ventriculomegalia leve preliminar";
    elseif d <= 15
        diagnostico = "Ventriculomegalia moderada preliminar";
    else
        diagnostico = "Ventriculomegalia severa preliminar";
    end

end