function history = fmincon_HE_2D
tic
%% SET UP SHARED VARIABLES WITH OUTFUN
history.x = [];
history.fval = [];
%% MATERIAL PROPERTIES
nely = 40; 
nelx = 40; 
volfrac = 0.3;
penal = 3; 
rmin = 1.4;
k0 = 1;           % Young's modulus of solid material
kmin = 1e-3;      % Young's modulus of void-like material
videocut = false;
%% SOLUTION INITIALIZATION
% VOLUME FRACTION
x0(1:nelx,1:nely) = volfrac;
xm = x0(:);
eq = ones(1,nelx*nely);
totalvolfrac = nely*nelx*volfrac;
% UPPER AND LOWER BOUNDS 
nele = nelx*nely;
    lb(1:nele) = 0.0;
    ub(1:nele) = 1.0;
%% DEFINE HEAT GENERATION AND BOUNDARY CONDITION
il = int32(nelx/2-nelx/20:nelx/2+nelx/20); jl = nely;    % Coordinates
fixedid = il*(nely+1)+(nely+1-jl);                       % Coordinates
fixeddof = reshape(fixedid,[],1);
nele = nelx*nely;
ndof = (nelx+1)*(nely+1);
F = sparse(1:ndof,1,0.01,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk;
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
edofVec = nodeids(:)+1;
edofMat = repmat(edofVec,1,4)+ repmat([0 nely+1 nely -1],nele,1);
iK = reshape(kron(edofMat,ones(4,1))',4*4*nele,1);
jK = reshape(kron(edofMat,ones(1,4))',4*4*nele,1);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% CALLBACK OBJECTIVE FUNCTION AND GRADIENT OBJECTIVE
objgrad = @(x)funcfdf2D(nelx,nely,penal,k0,kmin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,x);
hessfunc = @(x,lamda)hessinterior(nely,nelx,penal,k0,kmin,edofMat,KE,iK,jK,freedofs,U,F,nele,x);
options = optimoptions('fmincon','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,...
    'SubproblemAlgorithm','cg','TolFun',1e-6,'TolCon',1e-8,'MaxIteration',1000,...
    'MaxFunctionEvaluations', 10000,'HessianFcn',hessfunc)
% Run FMINCON.
delete FMINCON_2D_Data.txt
[x,fval] = fmincon(objgrad,xm,[],[],eq,totalvolfrac,lb,ub,[],options);
FMINCONtime = toc;
save('FMINCON_2D_Time.mat','FMINCONtime');
% HISTORICAL COMPLIANCE DATA
function stop = outfun(x,optimValues,state)
     stop = false;
     switch state
         case 'iter'
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x x];
     end
end
%% PLOT DENSITIES
iter = size(history.x);
 for i = 1:iter(2)-1
     x2 = vec2mat(history.x(:,i)',nelx);
     xmat = 1-x2;
     change = sum(abs(history.x(:,i+1)-history.x(:,i)));
     if change < 0.01
        videocut = true;
     end
     if videocut == false
        figure(1); imagesc(x2),colormap(flipud(gray)),caxis([0 1]);axis equal; axis off; P(i)=getframe(gcf);drawnow;
     end
 end
%% CREATE VIDEO WRITER & SAVE TOPOLOGY
 figure(2); imagesc(x2),colormap(flipud(gray)),caxis([0 1]);axis equal; axis off; drawnow;
 HE2DVid = VideoWriter('FMINCON_HE_2D.avi');
 HE2DVid.FrameRate = 20;
 open(HE2DVid);
 for i=1:length(P)
 frame = P(i) ;
 writeVideo(HE2DVid, frame);
 end
 close(HE2DVid);
 save('FMINCON_2D_Topology.mat','xmat');
end
% ----------------------------------------------------------------------
function [c,dcvect] = funcfdf2D(nelx,nely,penal,k0,kmin,U,edofMat,KE,iK,jK,freedofs,F,H,Hs,x)
warning('off','all');
xPhys = vec2mat(x,nelx);
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(kmin+xPhys(:)'.^penal*(k0-kmin)),16*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((kmin+xPhys.^penal*(k0-kmin)).*ce));
  dc = -penal*(k0-kmin)*xPhys.^(penal-1).*ce;
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  dc(:) = H*(xPhys(:).*dc(:))./Hs./max(1e-3,xPhys(:));
  dctrans = dc';
  dcvect = dctrans(:);
%% SAVING COMPLIANCE DATA AND CPU TIME
FEtime = toc;
diary FMINCON_2D_Data.txt
fprintf("%f",c)
fprintf("\n")
diary off 
end
function Hess = hessinterior (nely,nelx,penal,k0,kmin,edofMat,KE,iK,jK,freedofs,U,F,nvar,x)
xPhys = vec2mat(x,nelx);
% FE-ANALYSIS
sK = reshape(KE(:)*(kmin+xPhys(:)'.^penal*(k0-kmin)),16*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
% SECOND DERIVATIVE
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
ddc = 2*(penal*(k0-kmin)*xPhys.^(penal-1)).^2.*(1./(kmin+xPhys.^penal*(k0-kmin))).*ce;

Hvect = zeros(nvar,1); a = zeros(nvar,1); m = 1;
for i = 1:nely
    for j = 1:nelx
            a(m) = m;
            Hvect(m) = ddc(i,j);
            m = m + 1;
        end
end
Hess = sparse(a,a,Hvect,nvar,nvar);
end
%% THERMAL CONDUCTIVITY MATRIX
function [KE]=lk
KE = [ 2/3 -1/6 -1/3 -1/6
       -1/6 2/3 -1/6 -1/3
       -1/3 -1/6 2/3 -1/6
       -1/6 -1/3 -1/6 2/3]; 
end