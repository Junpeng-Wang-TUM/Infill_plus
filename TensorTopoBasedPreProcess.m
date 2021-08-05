function preEmbeddedElements = TensorTopoBasedPreProcess(U,nelx,nely,edofMat,E0,nu,thickness)
	preEmbeddedElements = [];
	
	%%1. Preparation for topology analysis
	%%1.1 Initialize Mesh Info
	InitializeMeshInfo(nelx,nely,edofMat);

	%%1.2 Compute Cartesian Stress Field
	ComputeCartesianStress(U,E0,nu);
	
	%%2. Topology Analysis
	
	
	%%3. Extract the elements to be pre-embedded
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
	deShapeFuncs_ = DeShapeFunction(gaussIPs);
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

function dN = DeShapeFunction(paraCoords)	
	%% 	paraCoords = [
	%%		s1 s2 s3 ...
	%%		t1 t2 t3 ...
	%%		p1 p2 p3 ...
	%% ]
	s = paraCoords(:,1);
	t = paraCoords(:,2);
	dN1ds = -(1-t); dN2ds = 1-t; 	dN3ds = 1+t; dN4ds = -(1+t);
	dN1dt = -(1-s); dN2dt = -(1+s); dN3dt = 1+s; dN4dt = 1-s;
	
	dN = zeros(2*length(s), 4);
	dN(1:2:end,:) = 0.25*[dN1ds dN2ds dN3ds dN4ds];
	dN(2:2:end,:) = 0.25*[dN1dt dN2dt dN3dt dN4dt];	
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
