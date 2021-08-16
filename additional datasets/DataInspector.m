%% Inspect the Data sets associated with the paper

% Author: Junpeng Wang (junpeng.wang@tum.de)
% Date: 2021-08-05

clear all
clc

fileName = 'Ex2.txt';

%%1. Read data 
fid = fopen(fileName, 'r');
tmp = fscanf(fid, '%s %s %s', 3);
tmp = fscanf(fid, '%s', 1);
tmp = fscanf(fid, '%d %d', 3);
nelx_ = tmp(1); nely_ = tmp(2);
boundingBox_ = [0 0; nelx_ nely_];
tmp = fscanf(fid, '%s %s', 2);
numEles = fscanf(fid, '%d', 1);
eleList = fscanf(fid, '%d', [1 numEles])';
tmp = fscanf(fid, '%s %s', 2);
numFixedNodes = fscanf(fid, '%d', 1);
if numFixedNodes>0
	fixingCond_ = fscanf(fid, '%d %d %d', [3 numFixedNodes])';
end
tmp = fscanf(fid, '%s %s', 2);
numLoadedNodes = fscanf(fid, '%d', 1);	
if numLoadedNodes>0
	loadingCond_ = fscanf(fid, '%d %e %e', [3 numLoadedNodes])';
end
fclose(fid);

%%2. Produce Cartesian mesh
allEles = zeros(nelx_*nely_,1); allEles(eleList) = 1;
carEleMapBack_ = find(1==allEles);
numEles_ = length(carEleMapBack_);
carEleMapForward_ = zeros(nelx_*nely_,1);	
carEleMapForward_(carEleMapBack_) = (1:numEles_)';
nodenrs = reshape(1:(nelx_+1)*(nely_+1), 1+nely_, 1+nelx_); 
eNodVec = reshape(nodenrs(1:end-1,1:end-1)+1, nelx_*nely_, 1);
eNodMat_ = repmat(eNodVec(carEleMapBack_),1,4);
tmp = [0 nely_+[1 0] -1]; 
for ii=1:4
	eNodMat_(:,ii) = eNodMat_(:,ii) + repmat(tmp(ii), numEles_,1);
end	
carNodMapBack_ = unique(eNodMat_);
numNodes_ = length(carNodMapBack_);
numDOFs_ = 2*numNodes_;
carNodMapForward_ = zeros((nelx_+1)*(nely_+1),1);
carNodMapForward_(carNodMapBack_) = (1:numNodes_)';		
for ii=1:4
	eNodMat_(:,ii) = carNodMapForward_(eNodMat_(:,ii));
end
nodeCoords_ = zeros(numNodes_,2);
xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nelx_:boundingBox_(2,1);
ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/nely_:boundingBox_(1,2);		
tmp = reshape(repmat(xSeed, nely_+1, 1), (nelx_+1)*(nely_+1), 1);
nodeCoords_(:,1) = tmp(carNodMapBack_);
tmp = repmat(ySeed, 1, nelx_+1)';
nodeCoords_(:,2) = tmp(carNodMapBack_);

if numFixedNodes>0
	fixingCond_(:,1) = carNodMapForward_(fixingCond_(:,1));
end
if numLoadedNodes>0
	loadingCond_(:,1) = carNodMapForward_(loadingCond_(:,1));
end

%%3. Show Data
xPatchs = nodeCoords_(:,1); xPatchs = xPatchs(eNodMat_');
yPatchs = nodeCoords_(:,2); yPatchs = yPatchs(eNodMat_');
cPatchs = zeros(size(yPatchs));
hdMesh = patch(xPatchs, yPatchs, cPatchs); hold on;
set(hdMesh, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1, 'EdgeColor', 'None');

if numFixedNodes>0
	tarNodeCoord = nodeCoords_(fixingCond_(:,1),:);
	hd1 = plot(tarNodeCoord(:,1), tarNodeCoord(:,2), 'x', 'color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 10); hold on; 
end

if numLoadedNodes>0
	coordLoadedNodes = nodeCoords_(loadingCond_(:,1),:);
	amplitudesF = mean(boundingBox_(2,:)-boundingBox_(1,:))/5 * loadingCond_(:,2:3)./vecnorm(loadingCond_(:,2:3), 2, 2);
	hold on; hd2 = quiver(coordLoadedNodes(:,1), coordLoadedNodes(:,2), amplitudesF(:,1), ...
		amplitudesF(:,2), 0, 'Color', [1 0 0], 'LineWidth', 2, 'MaxHeadSize', 1, 'MaxHeadSize', 1); 		
end	
axis equal; axis tight; axis off

