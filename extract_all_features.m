% December 10 2018 - Fixed bug where objects too small to extract morph
% feats from would be included as zeros for all morph features
function [feats] = extract_all_features(bounds, mask)
%% Extracts features from nuclei/glandular bounds
% bounds is a struct containing the centroid and boundary info
% img is a color or grayscale image used for computed for haralick features
% img can be omitted if not using haralick featuers
% select can allow you to control which group of features to run

    fprintf('\nExtracting Graph Features...')
    tic
    feats{1} = extract_graph_feats(bounds);
    toc
    
    
    fprintf('\nExtracting Morph Features...')
    tic
    feats{2} = extract_morph_feats(bounds);
    toc
    
    fprintf('\nExtracting CGT Features...')
    tic
    feats{3} = extract_CGT_features_v2(bounds);
    toc


    fprintf('\nExtracting Cluster Graph Features...')
    tic
    feats{4} = extract_cluster_graph_feats(bounds);
    toc
    
    fprintf('\nExtracting Basic Shape Features...')
    tic
    feats{5}=getNucleiFeatures(mask);
    toc
    


function [graphfeats] = extract_graph_feats(bounds)
gb_r = [bounds.centroid_r];
gb_c = [bounds.centroid_c];

[graphfeats] = get_graph_features([gb_r]',[gb_c]');


function [morphfeats]=extract_morph_feats(bounds)
%% Morph
gb_r = {bounds.r};
gb_c = {bounds.c};

badPTC = [];
feats = zeros(length(gb_r),25);
for j = 1:length(gb_r)
    %     try
    if(numel([gb_r{j}]') > 4)
        if size(gb_r{1},1) == 1
            [feat] = morph_features([gb_r{j}]',[gb_c{j}]');
        else
            [feat] = morph_features([gb_r{j}],[gb_c{j}]);
        end
        
        feats(j,:) = feat;
    else
        badPTC = [badPTC, j];
    end

end

feats(badPTC,:) = []; %remove bad PTCs

% morphfeats = [nanmean(feats) nanstd(feats) nanmedian(feats) min(feats)./max(feats)];
morphfeats = [nanmean(feats) nanstd(feats) nanmedian(feats) prctile(feats,5)./prctile(feats,95)];


function [CGTfeats] = extract_CGT_feats(bounds)
%% CGT
a=0.5;
r=0.2;
[CGTfeats, c_matrix, info, feats, network, edges] = extract_CGT_features(bounds,a,r);

function [CCGfeats,feature_list] = extract_cluster_graph_feats(bounds)
%CCG
info.alpha = 0.5;
info.radius = 0.2;
% build graph
alpha = info.alpha;
r = info.radius;
[VX,VY,x,y,edges] = construct_ccgs_optimized(bounds,alpha, r);
[CCGfeats,feature_list] = cluster_graph_features_optimized(bounds, edges);

function [Texturefeats haralickImages] = extract_texture_feats(img,mask)

%% texture 

%texture features are NOT used for PTC experiments

% alternatively 13 haralick
%addpath(genpath('Z:\Datasets\SatishData\haralick'))

info.dist = 1;
info.win = 1;
info.grays = 256;

n = 0;

if ndims(img) < 3
    gimg = img;
elseif ndims(img) == 3
    gimg = rgb2gray(img); % assume img is rgb
else
    fprintf('Unidentified image format')
end

if(isempty(mask))
    mask(:,:,1) = (gimg ~= max(max(gimg)));
end
%for grays = [64,128,256]
%    for win = [1,2,3]
%        n = n + 1
grays = info.grays;
win = info.win;
dist = info.dist;

himg = uint16(rescale_range(gimg,0,grays-1));
%himg = rescale_range(gimg,0,grays-1);
%f = haralick_img(himg,mask,grays,win,dist,1);
f = haralick_img(himg,mask,grays,win,dist,1);
Hf = f.img3;
%        HaralickFeat{n} = single(Hf);
%    end
%end

haralickImages = Hf;

for i = 1:size(Hf,3)
    feat = Hf(:,:,i);
    img_mask = mask(:,:,1);
    roi = feat(img_mask ==1);
    
    Texturefeats(n+1) = mean(roi(:));
    Texturefeats(n+2) = std(roi(:));
    %Texturefeats(n+3) = mode(roi(:));
    n = n + 2;
end


% count = 1;
% modifier = [{'mean '} {'standard deviation '}]
% for j = 1:numel(feat.names)
% for i = 1:numel(modifier)
%     HaralickTextureFeatureDescription{count} = [modifier{i} 'intensity ' feat.names{j}];
%     count = count + 1;
% end
% end
function [nucleiFeatures]=getNucleiFeatures(mask)

mask=logical(mask(:,:,1));

regionProperties = regionprops(mask,'Centroid', 'Area', 'Perimeter', 'Eccentricity', 'MinorAxisLength', ...
    'MajorAxisLength', 'Extent', 'Orientation','BoundingBox');

nucleiCentroids = cat(1, regionProperties.Centroid);

nucleiNum=size(regionProperties,1);
% for i=1:nucleiNum
%     nucleus=regionProperties(i);
%     bbox = nucleus.BoundingBox;
%     bbox = [round(bbox(1)) round(bbox(2)) (bbox(3) - 1) (bbox(4) - 1)];
%     roi = image(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3), :);
%     %groi = rgb2gray(roi);
% 
%     
% end

longness=[regionProperties.MajorAxisLength]./[regionProperties.MinorAxisLength];

Features = horzcat([regionProperties.Area]',[regionProperties.Eccentricity]',[regionProperties.MajorAxisLength]', ...
   [regionProperties.MinorAxisLength]', longness', [regionProperties.Extent]', [regionProperties.Orientation]'...
);

labels={'Area','Eccentricity', ...
    'Major Axis', 'Minor Axis','Major Axis / Minor Axis', 'Extent', 'Orientation'};
nucleiFeatures = [nanmean(Features) nanstd(Features)];


