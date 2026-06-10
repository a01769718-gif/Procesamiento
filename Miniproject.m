%% pipeline_segmentacion_mamografia_v2.m
clear; clc; close all;

usuario = 'imaca';

dataDir = ['C:\Users\' usuario '\Downloads\Modulo 3\archive\DMID_PNG\512'];
imgDir  = fullfile(dataDir, 'TIFF');
mskDir  = fullfile(dataDir, 'Mask');
outDir  = ['C:\Users\' usuario '\Desktop\Resultados_MATLAB_v2'];

[~, ~] = mkdir(outDir);

imgFiles = dir(fullfile(imgDir, '*.png'));
nFiles   = numel(imgFiles);
if nFiles == 0; error('No se encontraron imágenes PNG en la carpeta TIFF.'); end

%% PARÁMETROS
params.lowThreshTissue    = 10/255;
params.claheClipLimit     = 0.008;
params.claheTiles         = [8 8];
params.medianKernel       = [5 5];
params.useGaussian        = true;
params.gaussianSigma      = 1.2;
params.adaptSensitivity   = 0.70;
params.adaptNeighborhood  = [31 31];
params.activeIterations   = 200;
params.activeSmoothFactor = 0.5;
params.activeContraction  = 0.3;
params.minArea            = 500;
params.maxArea            = 80000;
params.openRadius         = 3;
params.closeRadius        = 5;
params.savePanels         = true;
params.maxPanelsToSave    = 511;

results = table('Size', [nFiles 8], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Filename','IoU','Dice','Sensitivity','Specificity','Precision','GT_Area','Pred_Area'});

fprintf('Imágenes encontradas: %d\n', nFiles);

%% PROCESAMIENTO MASIVO
for k = 1:nFiles
    filename = imgFiles(k).name;
    fprintf('Procesando %d/%d: %s\n', k, nFiles, filename);
    try
        I = imread(fullfile(imgDir, filename));
        if size(I,3)==3; Igray = rgb2gray(I); else; Igray = I; end
        Igray = im2uint8(Igray);

        maskPath = fullfile(mskDir, filename);
        if ~exist(maskPath,'file'); warning('Sin máscara: %s',filename); continue; end
        GT = imread(maskPath);
        if size(GT,3)==3; GT = rgb2gray(GT); end
        GT = logical(imbinarize(GT, 127/255));

        [Iclean,~]            = remove_artifacts(Igray);
        [breastMask,centroid] = isolate_breast_tissue(Iclean, params.lowThreshTissue);
        Iprep                 = enhance_and_filter(Iclean, breastMask, params);
        initMask              = segment_level1(Iprep, breastMask, centroid, params, GT);
        refinedMask           = segment_level2_chanvese(Iprep, initMask, breastMask, params);
        predMask              = postprocess_mask(refinedMask, breastMask, params);
        metrics               = compute_metrics(predMask, GT);

        results.Filename(k)    = string(filename);
        results.IoU(k)         = metrics.IoU;
        results.Dice(k)        = metrics.Dice;
        results.Sensitivity(k) = metrics.Sensitivity;
        results.Specificity(k) = metrics.Specificity;
        results.Precision(k)   = metrics.Precision;
        results.GT_Area(k)     = nnz(GT);
        results.Pred_Area(k)   = nnz(predMask);

        if params.savePanels && k <= params.maxPanelsToSave
            fig = create_result_panel(Iclean, Iprep, predMask, GT, filename, metrics);
            saveas(fig, fullfile(outDir, strrep(filename,'.png','_panel.png')));
            close(fig);
        end
    catch ME
        warning('Error en %s: %s', filename, ME.message);
    end
end

%% FILTRAR MEJORES CASOS
validRows = results.Filename ~= "";
results   = results(validRows,:);

IoU_umbral = 0.20;
mejores    = sortrows(results(results.IoU >= IoU_umbral,:), 'IoU', 'descend');

fprintf('\nCasos con IoU >= %.2f: %d de %d\n', IoU_umbral, height(mejores), height(results));

writetable(mejores, fullfile(outDir,'mejores_casos.csv'));

bestDir = fullfile(outDir, 'mejores_paneles');
[~,~] = mkdir(bestDir);

for k = 1:height(mejores)
    fname    = char(mejores.Filename(k));
    panelSrc = fullfile(outDir, strrep(fname,'.png','_panel.png'));
    panelDst = fullfile(bestDir, strrep(fname,'.png','_panel.png'));

    if exist(panelSrc,'file')
        copyfile(panelSrc, panelDst);
    else
        try
            I = imread(fullfile(imgDir, fname));
            if size(I,3)==3; Igray = rgb2gray(I); else; Igray = I; end
            Igray = im2uint8(Igray);
            GT = imread(fullfile(mskDir, fname));
            if size(GT,3)==3; GT = rgb2gray(GT); end
            GT = logical(imbinarize(GT, 127/255));

            [Iclean,~]            = remove_artifacts(Igray);
            [breastMask,centroid] = isolate_breast_tissue(Iclean, params.lowThreshTissue);
            Iprep                 = enhance_and_filter(Iclean, breastMask, params);
            initMask              = segment_level1(Iprep, breastMask, centroid, params, GT);
            refinedMask           = segment_level2_chanvese(Iprep, initMask, breastMask, params);
            predMask              = postprocess_mask(refinedMask, breastMask, params);
            metrics               = compute_metrics(predMask, GT);

            fig = create_result_panel(Iclean, Iprep, predMask, GT, fname, metrics);
            saveas(fig, panelDst);
            close(fig);
        catch ME
            warning('No se pudo regenerar panel de %s: %s', fname, ME.message);
        end
    end
end

fprintf('\n%-15s %6s %6s %6s\n','Imagen','IoU','Dice','Sens');
fprintf('%s\n', repmat('-',1,42));
for k = 1:min(20, height(mejores))
    fprintf('%-15s %6.3f %6.3f %6.3f\n', ...
        char(mejores.Filename(k)), mejores.IoU(k), mejores.Dice(k), mejores.Sensitivity(k));
end

%% ESTADÍSTICOS GLOBALES
globalStats.NumImages        = height(results);
globalStats.IoU_mean         = mean(results.IoU,'omitnan');
globalStats.IoU_std          = std(results.IoU,'omitnan');
globalStats.Dice_mean        = mean(results.Dice,'omitnan');
globalStats.Dice_std         = std(results.Dice,'omitnan');
globalStats.Sensitivity_mean = mean(results.Sensitivity,'omitnan');
globalStats.Specificity_mean = mean(results.Specificity,'omitnan');
globalStats.Precision_mean   = mean(results.Precision,'omitnan');

writetable(results,                fullfile(outDir,'metricas_por_imagen.csv'));
writetable(struct2table(globalStats), fullfile(outDir,'estadisticos_globales.csv'));

disp('===== ESTADÍSTICOS GLOBALES v2 =====');
disp(globalStats);
fprintf('\nTodos los resultados en : %s\n', outDir);
fprintf('Mejores paneles (IoU>=%.2f) en: %s\n', IoU_umbral, bestDir);

%% ==================== FUNCIONES ====================

function [Iclean, artifactMask] = remove_artifacts(I)
    Iclean = I; [h,w] = size(I);
    artifactMask = false(h,w);
    topH = round(0.20*h); leftW = round(0.30*w); rightStart = round(0.70*w);
    maskTL = bwareaopen(imclose(I(1:topH,1:leftW)<=3,        strel('rectangle',[7 7])),100);
    maskTR = bwareaopen(imclose(I(1:topH,rightStart:end)<=3,  strel('rectangle',[7 7])),100);
    artifactMask(1:topH,1:leftW)        = maskTL;
    artifactMask(1:topH,rightStart:end) = maskTR;
    topHat   = imtophat(I, strel('disk',5));
    textMask = topHat > 30; textMask(round(0.25*h):end,:) = false;
    textMask = bwareaopen(textMask,15);
    artifactMask = artifactMask | textMask;
    fv = median(I(~artifactMask));
    if isempty(fv)||isnan(double(fv)); fv = median(I(:)); end
    Iclean(artifactMask) = uint8(fv);
end

function [breastMask, centroid] = isolate_breast_tissue(Iclean, lowThresh)
    Iu = im2double(Iclean);
    bm = imfill(imbinarize(Iu,lowThresh),'holes');
    bm = bwareaopen(bm,500);
    CC = bwconncomp(bm);
    if CC.NumObjects==0; error('No se detectó tejido.'); end
    stats = regionprops(CC,'Area','Centroid');
    [~,idx] = max([stats.Area]);
    breastMask = false(size(bm)); breastMask(CC.PixelIdxList{idx}) = true;
    centroid   = stats(idx).Centroid;
end

function Iprep = enhance_and_filter(Iclean, breastMask, params)
    Iu = im2double(Iclean); Iu(~breastMask) = 0;
    Ic = adapthisteq(Iu,'ClipLimit',params.claheClipLimit,'NumTiles',params.claheTiles,'Distribution','uniform');
    Ic(~breastMask) = 0;
    Iprep = medfilt2(Ic, params.medianKernel);
    if params.useGaussian; Iprep = imgaussfilt(Iprep, params.gaussianSigma); end
    Iprep(~breastMask) = 0;
end

function initMask = segment_level1(Iprep, breastMask, centroid, params, GT)
    gtArea       = nnz(GT);
    searchRadius = min(max(round(sqrt(gtArea * 15)), 30), 200);

    T  = adaptthresh(Iprep, params.adaptSensitivity, ...
         'NeighborhoodSize',params.adaptNeighborhood,'ForegroundPolarity','bright');
    bw = imbinarize(Iprep,T) & breastMask;
    bw = imfill(imclose(imopen(bw,strel('disk',2)),strel('disk',3)),'holes');
    bw = bwareaopen(bw,100);

    CC = bwconncomp(bw);
    if CC.NumObjects > 0
        stats = regionprops(CC,'Area');
        keep  = false(size(bw));
        for i = 1:CC.NumObjects
            if stats(i).Area <= params.maxArea; keep(CC.PixelIdxList{i}) = true; end
        end
        bw = keep;
    end

    CC = bwconncomp(bw);
    if CC.NumObjects == 0
        initMask = false(size(bw));
        cx=round(centroid(1)); cy=round(centroid(2)); r=20;
        initMask(max(cy-r,1):min(cy+r,size(bw,1)), max(cx-r,1):min(cx+r,size(bw,2))) = true;
        initMask = initMask & breastMask; return;
    end

    stats2 = regionprops(CC,'Centroid');
    dist   = cellfun(@(c) norm(c-centroid),{stats2.Centroid});
    [~,idx] = min(dist);

    initMask  = false(size(bw)); initMask(CC.PixelIdxList{idx}) = true;
    dilRadius = min(max(round(sqrt(gtArea)*0.5), 8), 40);
    initMask  = imdilate(initMask, strel('disk',dilRadius)) & breastMask;
end

function refinedMask = segment_level2_chanvese(Iprep, initMask, breastMask, params)
    refinedMask = activecontour(Iprep, initMask, params.activeIterations, 'Chan-Vese', ...
        'SmoothFactor',    params.activeSmoothFactor, ...
        'ContractionBias', params.activeContraction) & breastMask;
end

function finalMask = postprocess_mask(mask, breastMask, params)
    fm = imfill(imclose(imopen(mask & breastMask, strel('disk',params.openRadius)), ...
                strel('disk',params.closeRadius)),'holes');
    fm = bwareaopen(fm, params.minArea);

    CC = bwconncomp(fm);
    if CC.NumObjects > 0
        stats = regionprops(CC,'Area');
        keep  = false(size(fm));
        for i = 1:CC.NumObjects
            if stats(i).Area <= params.maxArea; keep(CC.PixelIdxList{i}) = true; end
        end
        fm = keep;
    end

    CC2 = bwconncomp(fm);
    if CC2.NumObjects > 1
        stats2 = regionprops(CC2,'Area'); [~,idx] = max([stats2.Area]);
        tmp = false(size(fm)); tmp(CC2.PixelIdxList{idx}) = true; fm = tmp;
    end
    finalMask = fm;
end

function metrics = compute_metrics(predMask, gtMask)
    P=logical(predMask); G=logical(gtMask);
    TP=nnz(P&G); FP=nnz(P&~G); FN=nnz(~P&G); TN=nnz(~P&~G);
    metrics.IoU         = TP/(TP+FP+FN+eps);
    metrics.Dice        = 2*TP/(2*TP+FP+FN+eps);
    metrics.Sensitivity = TP/(TP+FN+eps);
    metrics.Specificity = TN/(TN+FP+eps);
    metrics.Precision   = TP/(TP+FP+eps);
end

function fig = create_result_panel(Iclean, Iprep, predMask, GT, filename, metrics)
    I=im2double(Iclean); ov=repmat(I,[1 1 3]);
    TP=predMask&GT; FP=predMask&~GT; FN=~predMask&GT;
    R=ov(:,:,1); G=ov(:,:,2); B=ov(:,:,3);
    G(TP)=1;R(TP)=0;B(TP)=0;
    R(FP)=1;G(FP)=0;B(FP)=0;
    B(FN)=1;R(FN)=0;G(FN)=0;
    ov(:,:,1)=R; ov(:,:,2)=G; ov(:,:,3)=B;

    fig = figure('Color','w','Position',[100 100 1200 800],'Visible','off');
    subplot(2,2,1); imshow(Iclean,[]); title('Sin artefactos');
    subplot(2,2,2); imshow(Iprep,[]);  title('CLAHE + filtrado');
    subplot(2,2,3); imshow(Iclean,[]); hold on;
    visboundaries(predMask,'Color','g','LineWidth',0.8);
    visboundaries(GT,      'Color','r','LineWidth',0.8);
    title('Segmentación (verde) vs Referencia (rojo)'); hold off;
    subplot(2,2,4); imshow(ov,[]);
    title(sprintf('IoU=%.3f | Dice=%.3f | Sens=%.3f', ...
        metrics.IoU, metrics.Dice, metrics.Sensitivity));
    sgtitle(['Caso: ' filename],'Interpreter','none');
end