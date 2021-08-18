function preEmbeddedElements = TensorTopoBasedPreProcess(U,nelx,nely,edofMat,E0,nu,thickness)
	preEmbeddedElements = [];
	
	%%1. Preparation for topology analysis
	%%1.1 Initialize Mesh Info
	InitializeMeshInfo(nelx,nely,edofMat);

	%%1.2 Compute Cartesian Stress Field
	ComputeCartesianStress(U,E0,nu);
	
	%%2. Topology Analysis
	TensorFieldTopologyAnalysis();
	
	%%3. Extract the elements to be pre-embedded
	preEmbeddedElements = ExtractElementsNearTopoSkelesTrisectorDegePoints(thickness);
	
	%%4. Show Topology of Stress Tensor Field
	figure(100); ShowStressTensorTopology();
end

function InitializeMeshInfo(nelx,nely,edofMat)
	global boundingBox_;
	global nelx_;
	global nely_;
	global eleSize_;
	global numEles_;
	global numNodes_;
	global nodeCoords_;
	global eNodMat_;
	global eDofMat_;
	global numNodsAroundEleVec_;
	global nodesOnBoundary_;
	global elementsOnBoundary_;

	global meshState_;	
	global carEleMapBack_; 
	global carEleMapForward_;
	global carNodMapBack_; 
	global carNodMapForward_;
	
	nelx_ = nelx;
	nely_ = nely;
	boundingBox_ = [0 0; nelx_ nely_];
	eleSize_ = 1;
	numEles_ = nelx_*nely_;
	numNodes_ = (nelx_+1)*(nely_+1);
	nodeCoords_ = zeros(numNodes_,2);

	xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nelx_:boundingBox_(2,1);
	ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/nely_:boundingBox_(1,2);		
	nodeCoords_(:,1) = reshape(repmat(xSeed, nely_+1, 1), (nelx_+1)*(nely_+1), 1);
	nodeCoords_(:,2) = repmat(ySeed, 1, nelx_+1)';	
	
	eDofMat_ = edofMat;
	eNodMat_ = eDofMat_(:,2:2:8)/2;
	
	numNodsAroundEleVec_ = zeros(numNodes_,1);
	for ii=1:numEles_
		iEleNodes = eNodMat_(ii,:);
		numNodsAroundEleVec_(iEleNodes,1) = numNodsAroundEleVec_(iEleNodes,1) + 1;
	end
	
	nodesOnBoundary_ = find(numNodsAroundEleVec_<4);
	allNodes = zeros(numNodes_,1); allNodes(nodesOnBoundary_) = 1;
	tmp = zeros(numEles_,1);
	for ii=1:4
		tmp = tmp + allNodes(eNodMat_(:,ii));
	end
	elementsOnBoundary_ = find(tmp>0);
	
	%%Only for the Specific 'Rectangle' Design Domain
	%%Further Work Needed when Considering Arbitrary Shapes
	meshState_ = ones(numEles_,1);
	carEleMapBack_ = (1:numEles_)';
	carEleMapForward_ = (1:numEles_)';
	carNodMapBack_ = (1:numNodes_)';
	carNodMapForward_ = (1:numNodes_)';
end

function ComputeCartesianStress(U, E0, nu)
	global numNodes_;
	global numEles_;
	global eNodMat_;
	global eDofMat_;
	global numNodsAroundEleVec_;
	global detJ_;
	global invJ_;	
	global deShapeFuncs_;
	global cartesianStressField_;
	
	%%1. Preparation
	gaussIPs = GaussianIntegral();
	deShapeFuncs_ = DeShapeFunction(gaussIPs(:,1), gaussIPs(:,2));
	[detJ_, invJ_] = CalcJacobi();
	iMatrixD = ElementElasticityMatrix(E0, nu);
	iMatrixB = ElementStrainMatrix(deShapeFuncs_, invJ_);	
	OTP = OuterInterpolationMat();
	
	%%2. Compute Cartesian Stress
	cartesianStressField_ = zeros(numNodes_, 3);
	for ii=1:numEles_
		iEleU = U(eDofMat_(ii,:),1);
		stressGaussPoints = iMatrixD * (iMatrixB*iEleU);
		stressNodes = OTP*stressGaussPoints;
		iNodes = eNodMat_(ii,:);
		cartesianStressField_(iNodes,:) = reshape(stressNodes,3,4)' + cartesianStressField_(iNodes,:);
	end
	cartesianStressField_ = cartesianStressField_./numNodsAroundEleVec_;
end

function [detJ, invJ] = CalcJacobi()
	%% only for 1st-order quad element
	global eNodMat_;
	global nodeCoords_;
	global deShapeFuncs_;
	nEND = 2;
	nEGIP = 4;
	detJ = zeros(nEGIP,1);
	invJ = sparse(nEND*nEGIP,nEND*nEGIP);	
	
	probeEleNods = nodeCoords_(eNodMat_(1,:)',:);
	for kk=1:nEGIP
		Jac = deShapeFuncs_(nEND*(kk-1)+1:nEND*kk,:)*probeEleNods;
		detJ(kk) = det(Jac);
		invJ(nEND*(kk-1)+1:nEND*kk, nEND*(kk-1)+1:nEND*kk) = inv(Jac);		
	end	
end

function N = ShapeFunction(s, t)
	%				   	   __s (parametric coordinate system)
	%				  	 /-t
	%				*4			*3
	%			*1			*2
	%
	%				nodes
	s = s(:);
	t = t(:);
	N = zeros(size(s,1), 4);
	N(:,1) = 0.25*(1-s).*(1-t);
	N(:,2) = 0.25*(1+s).*(1-t);
	N(:,3) = 0.25*(1+s).*(1+t);
	N(:,4) = 0.25*(1-s).*(1+t);	
end

function dN = DeShapeFunction(s, t)	
	s = s(:);
	t = t(:);
	dN1ds = -(1-t); dN2ds = 1-t; 	dN3ds = 1+t; dN4ds = -(1+t);
	dN1dt = -(1-s); dN2dt = -(1+s); dN3dt = 1+s; dN4dt = 1-s;
	
	dN = zeros(2*length(s), 4);
	dN(1:2:end,:) = 0.25*[dN1ds dN2ds dN3ds dN4ds];
	dN(2:2:end,:) = 0.25*[dN1dt dN2dt dN3dt dN4dt];	
end

function d2Shape = De2ShapeFunction(s, t)	
	s = s(:);
	t = t(:);
	numCoord = length(s);
	dN1dss = 0; dN2dss = 0; dN3dss = 0; dN4dss = 0;
	dN1dtt = 0; dN2dtt = 0; dN3dtt = 0; dN4dtt = 0;
	dN1dst = 0.25; dN2dst = -0.25; dN3dst = 0.25; dN4dst = -0.25;	
	
	d2Shape = repmat(s,3,4);
	d2Shape(1:3:end,:) = repmat([dN1dss	dN2dss	dN3dss	dN4dss], numCoord, 1);
	d2Shape(2:3:end,:) = repmat([dN1dtt	dN2dtt	dN3dtt	dN4dtt], numCoord, 1);
	d2Shape(3:3:end,:) = repmat([dN1dst	dN2dst	dN3dst	dN4dst], numCoord, 1);
end


function B = ElementStrainMatrix(dShape, invJ)
	derivatives = invJ * dShape;
	dNds1 = derivatives(1,:);	dNdt1 = derivatives(2,:);
	dNds2 = derivatives(3,:);	dNdt2 = derivatives(4,:);
	dNds3 = derivatives(5,:);	dNdt3 = derivatives(6,:);
	dNds4 = derivatives(7,:);	dNdt4 = derivatives(8,:);
	
	B11 = [dNds1(1) 0; 0 dNdt1(1); dNdt1(1) dNds1(1)];
	B12 = [dNds1(2) 0; 0 dNdt1(2); dNdt1(2) dNds1(2)];
	B13 = [dNds1(3) 0; 0 dNdt1(3); dNdt1(3) dNds1(3)];
	B14 = [dNds1(4) 0; 0 dNdt1(4); dNdt1(4) dNds1(4)];
	
	B21 = [dNds2(1) 0; 0 dNdt2(1); dNdt2(1) dNds2(1)];
	B22 = [dNds2(2) 0; 0 dNdt2(2); dNdt2(2) dNds2(2)];
	B23 = [dNds2(3) 0; 0 dNdt2(3); dNdt2(3) dNds2(3)];
	B24 = [dNds2(4) 0; 0 dNdt2(4); dNdt2(4) dNds2(4)];
	
	B31 = [dNds3(1) 0; 0 dNdt3(1); dNdt3(1) dNds3(1)];
	B32 = [dNds3(2) 0; 0 dNdt3(2); dNdt3(2) dNds3(2)];
	B33 = [dNds3(3) 0; 0 dNdt3(3); dNdt3(3) dNds3(3)];
	B34 = [dNds3(4) 0; 0 dNdt3(4); dNdt3(4) dNds3(4)];
	
	B41 = [dNds4(1) 0; 0 dNdt4(1); dNdt4(1) dNds4(1)];
	B42 = [dNds4(2) 0; 0 dNdt4(2); dNdt4(2) dNds4(2)];
	B43 = [dNds4(3) 0; 0 dNdt4(3); dNdt4(3) dNds4(3)];
	B44 = [dNds4(4) 0; 0 dNdt4(4); dNdt4(4) dNds4(4)];
	
	B = [
		B11 B12 B13 B14
		B21 B22 B23 B24
		B31 B32 B33 B34
		B41 B42 B43 B44
	];	
end

function D = ElementElasticityMatrix(E, nu)	
	D = zeros(12);
	HL = zeros(3);
	HL(1,1) = E/(1-nu^2); HL(1,2) = E*nu/(1-nu^2);		
	HL(2,1) = HL(1,2); HL(2,2) = HL(1,1);	
	HL(3,3) = E/2/(1+nu);
	
	for ii=1:4
		index = (ii-1)*3+1:ii*3;
		D(index,index) = HL;
	end		
	D = sparse(D);
end

function gaussIPs = GaussianIntegral();
	sqrt33 = sqrt(3)/3; 
	gaussIPs = [-1 1 1 -1; -1 -1 1 1]' * sqrt33;
end

function outerInterpolationMatrix = OuterInterpolationMat()
	gaussIPs = GaussianIntegral();
	s = gaussIPs(:,1);
	t = gaussIPs(:,2);
	N = ShapeFunction(s,t);
	sFM = sparse(12,12);
	ii = 3*(1:4);
	sFM(1,ii-2) = N(1,:); sFM(2,ii-1) = N(1,:); sFM(3,ii) = N(1,:);
	sFM(4,ii-2) = N(2,:); sFM(5,ii-1) = N(2,:);	sFM(6,ii) = N(2,:);	
	sFM(7,ii-2) = N(3,:); sFM(8,ii-1) = N(3,:); sFM(9,ii) = N(3,:);	
	sFM(10,ii-2) = N(4,:); sFM(11,ii-1) = N(4,:); sFM(12,ii) = N(4,:);	
	outerInterpolationMatrix = inv(sFM);
end


function TensorFieldTopologyAnalysis()
	global potentialDegenerateElements_;
	global degeneratePoints_;
	
	%%1. identify potential 'degenerate' elements
	ExtractPotentialDegenerateElements();

	%%2. identify degenerate points
	IdentifyingDegeneratePoints(potentialDegenerateElements_);	
	
	%%3. post-processing degenerate points
	degeneratePoints_ = PostProcessDegeneratePoints();
	
	%%4. get the topology skeletons
	ComputeTopologicalSkeletons2D();	
end

function ExtractPotentialDegenerateElements()
	global numEles_;
	global eNodMat_;	
	global cartesianStressField_;
	global potentialDegenerateElements_;
	potentialDegenerateElements_ = [];
	for ii=1:numEles_
		eleStress = cartesianStressField_(eNodMat_(ii,:)',:);
		opt = DegenrationMeasure(eleStress);
		if 1==opt, potentialDegenerateElements_(end+1,1) = ii; end			
	end
end

function opt = DegenrationMeasure(tar)
	discriminants = DiscriminantConstraintFuncs(tar);
	v1 = discriminants(:,1);
	
	bool1_1 = v1(1)>0 && v1(2)>0 && v1(3)>0 && v1(4)>0;
	bool1_2 = v1(1)<0 && v1(2)<0 && v1(3)<0 && v1(4)<0;
	bool1 = bool1_1 || bool1_2;
	
	v2 = discriminants(:,2);
	bool2_1 = v2(1)>0 && v2(2)>0 && v2(3)>0 && v2(4)>0;
	bool2_2 = v2(1)<0 && v2(2)<0 && v2(3)<0 && v2(4)<0;
	bool2 = bool2_1 || bool2_2;
	
	if bool1 || bool2, opt = 0; else, opt = 1; end	
end

function discriminants = DiscriminantConstraintFuncs(eleStress)
	discriminants = [eleStress(:,1)-eleStress(:,2), eleStress(:,3)];
end

function IdentifyingDegeneratePoints(potentialElements)	
	global consideredDegenerateElements_; 
	global numConsideredDegenerateElements_; 
	global eNodMat_;
	global cartesianStressField_;
	global paraCoordListDegeEles_;
	
	consideredDegenerateElements_ = potentialElements(:);
	numConsideredDegenerateElements_ = length(consideredDegenerateElements_);	
	paraCoordListDegeEles_ = zeros(numConsideredDegenerateElements_,2);
	for ii=1:numConsideredDegenerateElements_
		iEleStress = cartesianStressField_(eNodMat_(potentialElements(ii),:)',:); 
		v1 = DiscriminantConstraintFuncs(iEleStress);
		[paraCoord, ~, ~, ~] = NewtonIteration(v1, zeros(1,size(v1,2)));
		paraCoordListDegeEles_(ii,:) = paraCoord;
	end
end

function [paraCoordinates, res, opt, index] = NewtonIteration(vtxVec, target)
	%% solving a nonlinear system of equasions by Newton-Rhapson's method
	%%	f1(s,t) = tar1
	%%	f2(s,t) = tar2
	opt = 0;
	normTar = norm(target);
	errThreshold = 1.0e-10; RF = 100*errThreshold;	
	s = -0.0; t = -0.0; maxIts = 150;
	index = 0;
	
	for ii=1:maxIts
		index = index+1;
		c0 = ShapeFunction(s, t)';
		dShape = DeShapeFunction(s, t);
		dns = dShape(1,:)';
		dnt = dShape(2,:)';
		d2Shape = De2ShapeFunction(s, t);
		dnss = d2Shape(1,:)';
		dntt = d2Shape(2,:)';
		dnst = d2Shape(3,:)';
		
		q = vtxVec' * c0;
		dqs = vtxVec' * dns;
		dqt = vtxVec' * dnt;

		dfdv1 = [dqs';dqt'];
		b = dfdv1*(q-target');
		if 0==normTar
			res = norm(q-target');
		else
			res = norm(b);
		end
		if res < errThreshold, break; end			
		
		dfdss = vtxVec'*dnss;
		dfdtt = vtxVec'*dntt;
		dfdst = vtxVec'*dnst;
		A11 = dfdss' * (q-target') + norm(dqs)^2;
		A22 = dfdtt' * (q-target') + norm(dqt)^2;
		A12 = dfdst' * (q-target') + dqs'*dqt;
		A21 = A12;
		A = [A11 A12; A21 A22]; x = A\(-b);		
		s = s + x(1); t = t + x(2);	
	end
	if res <= errThreshold && abs(s)<=RF+1 && abs(t)<=RF+1
		opt = 1;
	end
	paraCoordinates = [s t];
end

function extractedDegeneratePoints = PostProcessDegeneratePoints()
	global consideredDegenerateElements_;
	global paraCoordListDegeEles_;
	global numConsideredDegenerateElements_;
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global thresholdDPE_;
	global eleSize_;
	
	extractedDegeneratePoints = DegeneratePointStruct();
	RF = 1.0e-6; %%relaxation factor
	phyCoordList = [];
	index = 0;
	for ii=1:numConsideredDegenerateElements_	
		paraCoord = paraCoordListDegeEles_(ii,:);
		if abs(paraCoord(1))>RF+1 || abs(paraCoord(2))>RF+1, continue; end
		iDegePot = DegeneratePointStruct();
		iDegePot.eleIndex = consideredDegenerateElements_(ii);
		iDegePot.paraCoord = paraCoord;
		tarEleNodeIndices = eNodMat_(iDegePot.eleIndex,:)';
		iNodeCoord = nodeCoords_(tarEleNodeIndices,:);
		shapeFuncs = ShapeFunction(iDegePot.paraCoord(1), iDegePot.paraCoord(2));
		iDegePot.phyCoord = shapeFuncs*iNodeCoord;
		phyCoordList(end+1,:) = iDegePot.phyCoord;	
		iEleStress = cartesianStressField_(tarEleNodeIndices,:);
		iDegePot.cartesianStress = shapeFuncs*iEleStress;
		ps = ComputePrincipalStress(iDegePot.cartesianStress);
		iDegePot.principalStress = ps([4 1]);
		directDegenerancyMetric = abs(iDegePot.principalStress(1)-iDegePot.principalStress(2)) / ...
				abs(iDegePot.principalStress(1)+iDegePot.principalStress(2));	
		if directDegenerancyMetric>thresholdDPE_, continue; end
		iDegePot.directDegenerancyExtentMetric = directDegenerancyMetric;
		index = index + 1;	
		extractedDegeneratePoints(index,1) = iDegePot;				
	end
	if 0 == extractedDegeneratePoints(1).eleIndex, extractedDegeneratePoints(1) = []; end %% There is no degenerate point
end

function val = DegeneratePointStruct()
	vec = struct('ith', 0, 'length', 0, 'vec', [], 'index', []);
	PSL = PrincipalStressLineStruct();
	global eleType_;
	val = struct(	...
		'eleIndex',							0,	...
		'paraCoord',						[],	...		
		'phyCoord',							[],	...
		'cartesianStress',					[], ...
		'principalStress',					[],	...
		'directDegenerancyExtentMetric', 	[], ...
		'tangentList',						[],	...
		'stress2phy',						[],	...
		'abcd',								[],	...
		'delta',							0,	...
		'majorSkeletons',					PSL,...
		'minorSkeletons',					PSL,...
		'separatrices',						vec	...
	);
end

function val = PrincipalStressLineStruct()
	val = struct(...
		'ith',						0, 	...
		'length',					0,	...
		'midPointPosition',			0,	...		
		'phyCoordList',				[], ...
		'eleIndexList',				[], ...
		'principalStressList',		[] ...
	);	
end

function principalStress = ComputePrincipalStress(cartesianStress)
	principalStress = zeros(size(cartesianStress,1), 1+2+1+2);
	iPS = zeros(1, 6);
	for ii=1:size(cartesianStress,1)
		iCartesianStress = cartesianStress(ii,:);
		A = iCartesianStress([1 3; 3 2]);
		[eigenVec, eigenVal] = eig(A);
		iPS([1 4]) = diag(eigenVal);
		iPS([2 3 5 6]) = reshape(eigenVec,1,4);
		principalStress(ii,:) = iPS;
	end		
end

function ComputeTopologicalSkeletons2D()
	global boundingBox_;
	global eleSize_;
	global nodeCoords_; 
	global eNodMat_;
	global cartesianStressField_;
	global degeneratePoints_;
	global elementsOnBoundary_;
	global invJ_;
	global tracingStepWidth_;
	
	%%1. get the derivatives of cartesian stresses at degenerate points with respect to the cartesian coordinates
	%%  compute rotational invariant
	for ii=1:length(degeneratePoints_)
		s = degeneratePoints_(ii).paraCoord(1); t = degeneratePoints_(ii).paraCoord(2);
		dShape = DeShapeFunction(s, t);
		eleCoords = nodeCoords_(eNodMat_(degeneratePoints_(ii).eleIndex,:)',:);
		eleNodeCartesionStresses = cartesianStressField_(eNodMat_(degeneratePoints_(ii).eleIndex,:)',:);
		dN2dPhyC = invJ_(1:2,1:2)*dShape;	
		degeneratePoints_(ii).stress2phy = [(dN2dPhyC(1,:)*eleNodeCartesionStresses )' (dN2dPhyC(2,:)*eleNodeCartesionStresses )' ];
		
		a = (degeneratePoints_(ii).stress2phy(1,1) - degeneratePoints_(ii).stress2phy(2,1))/2; 
		b = (degeneratePoints_(ii).stress2phy(1,2) - degeneratePoints_(ii).stress2phy(2,2))/2; 
		c = degeneratePoints_(ii).stress2phy(3,1); 
		d = degeneratePoints_(ii).stress2phy(3,2);
		degeneratePoints_(ii).abcd = [a b c d];
		degeneratePoints_(ii).delta = a*d - b*c; %% Negative -> Trisector; Positive -> Wedge
		
		rawRoots = roots([d (c+2*b) (2*a-d) -c]);
		degeneratePoints_(ii).tangentList = rawRoots(0==imag(rawRoots));		
	end
	
	%%1.1 Exclude Degenerate Points located in the Boundary Elements
	%% "An Experience-based Solution to Ease up the Numerical Instability Caused by the POTENTIAL  Jaggy Boundary of Cartesian Mesh"
	elementsWithDegeneratePoints = [degeneratePoints_.eleIndex];
	[~, boundaryElementIndicesWithDegeneratePoints] = intersect(elementsWithDegeneratePoints, elementsOnBoundary_);
	degeneratePoints_(boundaryElementIndicesWithDegeneratePoints) = [];
	
	%%2. get the topology skeletons
	tracingStepWidth_ = eleSize_;
	stopCond = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/tracingStepWidth_);	
	for ii=1:length(degeneratePoints_)
		degeneratePoints_(ii).majorSkeletons = repmat(degeneratePoints_(ii).majorSkeletons, length(degeneratePoints_(ii).tangentList), 1);
		degeneratePoints_(ii).minorSkeletons = repmat(degeneratePoints_(ii).minorSkeletons, length(degeneratePoints_(ii).tangentList), 1);
		for jj=1:length(degeneratePoints_(ii).tangentList)
			seed = [degeneratePoints_(ii).eleIndex degeneratePoints_(ii).paraCoord];
			iniDir = [1 degeneratePoints_(ii).tangentList(jj) ]; iniDir = iniDir/norm(iniDir);
			[maxPSL, minPSL] = GeneratePrincipalStressLines(seed, iniDir, stopCond);
			
			dis0 = degeneratePoints_(ii).phyCoord - maxPSL.phyCoordList(1,:); dis0 = norm(dis0);
			dis1 = degeneratePoints_(ii).phyCoord - maxPSL.phyCoordList(end,:); dis1 = norm(dis1);		
			if dis0 < dis1
				maxPSL.eleIndexList = maxPSL.eleIndexList(maxPSL.midPointPosition:maxPSL.length,:);
				maxPSL.phyCoordList = maxPSL.phyCoordList(maxPSL.midPointPosition:maxPSL.length,:);
				maxPSL.principalStressList = maxPSL.principalStressList(maxPSL.midPointPosition:maxPSL.length,:);
				maxPSL.midPointPosition = 1;
				maxPSL.length = length(maxPSL.eleIndexList);
			else
				maxPSL.eleIndexList = flip(maxPSL.eleIndexList(1:maxPSL.midPointPosition,:));
				maxPSL.phyCoordList = flip(maxPSL.phyCoordList(1:maxPSL.midPointPosition,:));
				maxPSL.principalStressList = flip(maxPSL.principalStressList(1:maxPSL.midPointPosition,:)); 
				maxPSL.midPointPosition = 1;
				maxPSL.length = length(maxPSL.eleIndexList);
			end
			degeneratePoints_(ii).majorSkeletons(jj) = maxPSL;
			
			dis0 = degeneratePoints_(ii).phyCoord - minPSL.phyCoordList(1,:); dis0 = norm(dis0);
			dis1 = degeneratePoints_(ii).phyCoord - minPSL.phyCoordList(end,:); dis1 = norm(dis1);
			if dis0 < dis1
				minPSL.eleIndexList = minPSL.eleIndexList(minPSL.midPointPosition:minPSL.length,:);
				minPSL.phyCoordList = minPSL.phyCoordList(minPSL.midPointPosition:minPSL.length,:);
				minPSL.principalStressList = minPSL.principalStressList(minPSL.midPointPosition:minPSL.length,:);
				minPSL.midPointPosition = 1; minPSL.length = length(minPSL.eleIndexList);
			else
				minPSL.eleIndexList = flip(minPSL.eleIndexList(1:minPSL.midPointPosition,:));
				minPSL.phyCoordList = flip(minPSL.phyCoordList(1:minPSL.midPointPosition,:));
				minPSL.principalStressList = flip(minPSL.principalStressList(1:minPSL.midPointPosition,:)); 
				minPSL.midPointPosition = 1; minPSL.length = length(minPSL.eleIndexList);
			end		
			degeneratePoints_(ii).minorSkeletons(jj) = minPSL;	
		end
	end
end

function [majorPSL, minorPSL] = GeneratePrincipalStressLines(initialSeed, iniDir, limiSteps)
	majorPSL = PrincipalStressLineStruct();
	minorPSL = PrincipalStressLineStruct();
	
	%%1. Spot the Starting Point
	[eleIndex, phyCoord, principalStress] = PreparingForTracing(initialSeed);
	%%2. Compute PSL(s)
	psDir = [5 6];
	majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps);
	psDir = [2 3];
	minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps);
end

function iPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps)
	global tracingStepWidth_;
	iPSL = PrincipalStressLineStruct();
	PSLphyCoordList = phyCoord;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = principalStress;
	%% tracing along first direction (v1)
	nextPoint = phyCoord + tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, iniDir, eleIndex, psDir, limiSteps);
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	%% tracing along second direction (-v1)			
	nextPoint = phyCoord - tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, -iniDir, eleIndex, psDir, limiSteps);
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		eleIndexList = flip(eleIndexList);
		principalStressList = flip(principalStressList);				
	end
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	iPSL.midPointPosition = size(phyCoordList,1)+1;	
	%%2.3 finish Tracing the current major PSL			
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.principalStressList = PSLprincipalStressList;				
end

function [eleIndex, phyCoord, principalStress] = PreparingForTracing(initialSeed)
	global nodeCoords_; global eNodMat_;
	global cartesianStressField_;
	eleIndex = 0;
	phyCoord = 0; 
	principalStress = 0;
	formatedSeed = initialSeed;
	eleIndex = formatedSeed(1,1);					
	NIdx = eNodMat_(eleIndex,:)';
	eleNodeCoords = nodeCoords_(NIdx,:);
	eleCartesianStress = cartesianStressField_(NIdx,:);
	paraCoord = formatedSeed(1, 2:3);
	shapeFuncs = ShapeFunction(paraCoord(1), paraCoord(2));	
	phyCoord = shapeFuncs*eleNodeCoords;							
	interpolatedCartesianStress = shapeFuncs*eleCartesianStress;
	principalStress = ComputePrincipalStress(interpolatedCartesianStress);	
end

function [phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, iniDir, elementIndex, typePSL, limiSteps)
	%% Tracing the PSL by 2-nd order Runge-Kutta Scheme 
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global tracingStepWidth_; 
	
	phyCoordList = zeros(limiSteps,2);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,6);
	index = 0;	
	
	intgerScheme = 'RK2'; %% 'RK2', 'EULER'
	switch intgerScheme		
		case 'EULER'
			[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
			while 1==bool1
				index = index + 1; if index > limiSteps, index = index-1; break; end
				cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
				shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
				cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);					
				nextDir = DirectionSelecting(iniDir, principalStress(typePSL), -principalStress(typePSL));
					
				if 0 == AngleTerminationCondition(iniDir, nextDir), index = index-1; break; end	
				iniDir = nextDir;
				phyCoordList(index,:) = nextPoint;
				eleIndexList(index,:) = elementIndex;
				principalStressList(index,:) = principalStress;					
				nextPoint = nextPoint + tracingStepWidth_*iniDir;
				[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
			end				
		case 'RK2'
			%%initialize initial k1 and k2
			k1 = iniDir;
			iniPot = nextPoint - k1*tracingStepWidth_;
			midPot = nextPoint - k1*tracingStepWidth_/2;
			
				
			[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(midPot);
			if bool1
				cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
				shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
				cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
				k2 = DirectionSelecting(k1, principalStress(typePSL), -principalStress(typePSL));
				nextPoint = iniPot + tracingStepWidth_*k2;
				[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
				while 1==bool1
					index = index + 1; if index > limiSteps, index = index-1; break; end
					%%k1
					cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
					shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
					cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
					principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);					
					k1 = DirectionSelecting(iniDir, principalStress(typePSL), -principalStress(typePSL));	
					if 0 == AngleTerminationCondition(iniDir, k1), index = index-1; break; end
					%%k2
					midPot = nextPoint + k1*tracingStepWidth_/2;
					[elementIndex2, paraCoordinates2, bool1] = FindAdjacentElement(midPot);
					if ~bool1, index = index-1; break; end
					cartesianStress2 = cartesianStressField_(eNodMat_(elementIndex2,:)', :);
					shapeFuncs = ShapeFunction(paraCoordinates2(1), paraCoordinates2(2));
					cartesianStressOnGivenPoint2 = shapeFuncs*cartesianStress2;
					principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
					k2 = DirectionSelecting(k1, principalStress2(typePSL), -principalStress2(typePSL));		
					%%store	
					iniDir = k1;
					phyCoordList(index,:) = nextPoint;
					eleIndexList(index,:) = elementIndex;
					principalStressList(index,:) = principalStress;	
					%%next point
					nextPoint = nextPoint + tracingStepWidth_*k2;		
					[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
				end		
			end		
	end
	phyCoordList = phyCoordList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	principalStressList = principalStressList(1:index,:);				
end

function val = AngleTerminationCondition(dirct1, dirct2)
	angle = acos((dirct1*dirct2') / (norm(dirct1)*norm(dirct2)));
	if angle > pi/20
		val = 0;
	else
		val = 1;
	end
end

function targetDirection = DirectionSelecting(originalVec, Vec1, Vec2)
	normOriVec = norm(originalVec); normVec1 = norm(Vec1); normVec2 = norm(Vec2);
	angle1 = acos(originalVec*Vec1');
	angle2 = acos(originalVec*Vec2');
	if angle1 < angle2
		targetDirection = Vec1;
	else
		targetDirection = Vec2;
	end
end

function [nextElementIndex, paraCoordinates, opt] = FindAdjacentElement(physicalCoordinates)
	global nelx_; 
	global nely_; 
	global eleSize_;
	global nodeCoords_; 
	global eNodMat_;
	global meshState_; 
	global carEleMapForward_;
	global boundingBox_;
	Lbound = boundingBox_(1,:);
	nextElementIndex = 0; paraCoordinates = []; opt = 0;
	
	physicalCoordinates = physicalCoordinates - Lbound;
	if 0==physicalCoordinates(1)
		eleX = 1;				
	else
		eleX = ceil(physicalCoordinates(1)/eleSize_);
		if eleX<1 || eleX>nelx_, return; end
	end
	if 0==physicalCoordinates(2)
		eleY = 1;
	else
		eleY = ceil(physicalCoordinates(2)/eleSize_);
		if eleY<1 || eleY>nely_, return; end
	end				
	
	tarEle = nely_*(eleX-1)+(nely_-eleY+1);
	if meshState_(tarEle)	
		nextElementIndex = carEleMapForward_(tarEle);
		opt = 1;
		relatedNodes = eNodMat_(nextElementIndex,:);
		relatedNodeCoords = nodeCoords_(relatedNodes',:)-Lbound;
		paraCoordinates = 2*(physicalCoordinates - relatedNodeCoords(1,:)) / eleSize_ - 1;
	end
end

function preEmbeddedElements = ExtractElementsNearTopoSkelesTrisectorDegePoints(thickness)
	global degeneratePoints_;
	
	PSList = PrincipalStressLineStruct();
	
	%%1. Exclude the Wedge Degenerate points
	for ii=1:length(degeneratePoints_)
		if degeneratePoints_(ii).delta < 0
			PSList(end+1:end+6,1) = [degeneratePoints_(ii).majorSkeletons; degeneratePoints_(ii).minorSkeletons];
		end
	end
	PSList(1) = [];
	
	%%2. Compact PSLs
	tarIndice = [];
	for ii=1:length(PSList)
		PSList(ii).eleIndexList = PSList(ii).eleIndexList';
		if PSList(ii).length > 5 %%(Exclude the 'SHORT' PSLs)
			tarIndice(end+1,1) = ii;
		end
	end
	PSList = PSList(tarIndice);
	
	%%3. Extraction
	if isempty(PSList)
		preEmbeddedElements = [];
	else	
		preEmbeddedElements = [PSList.eleIndexList];
		index = 2;
		while index<=thickness
			preEmbeddedElements = RelateAdjacentElements(preEmbeddedElements);
			index = index + 1;
		end				
	end
end

function oEleList = RelateAdjacentElements(iEleList)
	global nelx_; global nely_;
	global carEleMapForward_;
	global carEleMapBack_;
	
	iEleListMapBack = carEleMapBack_(iEleList);
	%%	1	4	7
	%%	2	5*	8
	%%	3	6	9
	%%[eleX, eleY] = NodalizeDesignDomain([nelx_-1 nely_-1], [1 1; nelx_ nely_]);
	nx = nelx_-1; ny = nely_-1;
	dd = [1 1; nelx_ nely_];
	xSeed = dd(1,1):(dd(2,1)-dd(1,1))/nx:dd(2,1);
	ySeed = dd(2,2):(dd(1,2)-dd(2,2))/ny:dd(1,2);		
	eleX = repmat(xSeed, ny+1, 1); eleX = reshape(eleX, (nx+1)*(ny+1), 1);
	eleY = repmat(ySeed, 1, nx+1)';	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	eleX = eleX(iEleListMapBack);
	eleY = eleY(iEleListMapBack);
	tmpX = [eleX-1 eleX-1 eleX-1  eleX eleX eleX  eleX+1 eleX+1 eleX+1]; tmpX = tmpX(:);
	tmpY = [eleY+1 eleY eleY-1  eleY+1 eleY eleY-1  eleY+1 eleY eleY-1]; tmpY = tmpY(:);
	xNegative = find(tmpX<1); xPositive = find(tmpX>nelx_);
	yNegative = find(tmpY<1); yPositive = find(tmpY>nely_);
	allInvalidEles = unique([xNegative; xPositive; yNegative; yPositive]);
	tmpX(allInvalidEles) = []; tmpY(allInvalidEles) = [];
	oEleListMapBack = nely_*(tmpX-1) + nely_-tmpY + 1;
	
	oEleList = carEleMapForward_(oEleListMapBack);
	oEleList(oEleList<1) = []; oEleList = unique(oEleList);
end

function ShowStressTensorTopology()
	global numEles_;
	global eNodMat_;
	global nodeCoords_;
	global nodesOnBoundary_;
	global degeneratePoints_;
	
	if length(degeneratePoints_) > 0
		%%2. draw degenerate points and topological skeletons
		for ii=1:length(degeneratePoints_)		
			% 	
			for jj=1:length(degeneratePoints_(ii).majorSkeletons)
				if degeneratePoints_(ii).majorSkeletons(jj).length > 5
					plot(degeneratePoints_(ii).majorSkeletons(jj).phyCoordList(:,1), degeneratePoints_(ii).majorSkeletons(jj).phyCoordList(:,2), ...
						'-', 'color', [1 0 0], 'LineWidth', 3); hold on; %%[252 141 98]/255
				end
			end
			for jj=1:length(degeneratePoints_(ii).minorSkeletons)
				if degeneratePoints_(ii).minorSkeletons(jj).length > 5
					plot(degeneratePoints_(ii).minorSkeletons(jj).phyCoordList(:,1), degeneratePoints_(ii).minorSkeletons(jj).phyCoordList(:,2), ...
						'-', 'color', [0 0 1], 'LineWidth', 3); hold on %%[102 194 165]/255
				end
			end
			if degeneratePoints_(ii).delta > 0
				plot(degeneratePoints_(ii).phyCoord(1), degeneratePoints_(ii).phyCoord(2), 'sk', 'LineWidth', 4, 'MarkerSize', 20); hold on
			else
				plot(degeneratePoints_(ii).phyCoord(1), degeneratePoints_(ii).phyCoord(2), 'ok', 'LineWidth', 4, 'MarkerSize', 20); hold on
			end	
			
		end
	end

	%%Show silhouette
	edgeIndices = eNodMat_(:, [1 2 2 1  2 3 3 2  3 4 4 3  4 1 1 4])';
	edgeIndices = reshape(edgeIndices(:), 4, 4*numEles_);	
	tmp = zeros(size(nodeCoords_,1),1); tmp(nodesOnBoundary_) = 1;
	tmp = tmp(edgeIndices); tmp = sum(tmp,1);
	boundaryEleEdges = edgeIndices(:,find(4==tmp));
	xPatchs = nodeCoords_(:,1); xPatchs = xPatchs(boundaryEleEdges);
	yPatchs = nodeCoords_(:,2); yPatchs = yPatchs(boundaryEleEdges);		
	cPatchs = zeros(size(yPatchs));
	hd = patch(xPatchs, yPatchs, cPatchs); hold on;
	set(hd, 'FaceColor', 'None', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
	
	axis equal; axis tight; axis off;
end