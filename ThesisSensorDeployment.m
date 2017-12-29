% Thesis sensor deployment

clear;
% close all;

T = 1e1; % maximum time
nsteps = 1e1; % time interval
volfrac = 0.3; % volume fraction

% Set parameters
nelx = 4; % horizontal number of elements (left to right)
nely = 4; % vertical number of elements (top to down)
penal = 3; % polynomial order to define density-young's modulus relationship
E0 = 1; % young's modulus at density=1
Emin = 1e-9; % young's modulus at density=0, keep this small
M0 = 1;
Mmin = 1e-9;
nu = 0.3; % poisson's ratio


test = 1000;
DReal = zeros(1,test);
AvgShortDist = zeros(1,test);
Rstore = zeros(nelx,test);


for j = 1:test
% load structure.mat; % load xPhys

nf = nelx; % number of forces
num_observer = nf; % number of observers

% Left size of the beam is fixed to the ground
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
p = size(freedofs,2);

% set observer
S = zeros(num_observer,p); 
% S(:,randperm(p,num_observer))=eye(num_observer);
%S(:,2*(1:nelx)*(nely+1)) = eye(nelx); % [ORIGNAL] put y-axis sensors at the bottom of the beam
% R = randperm(p/2,num_observer); %creates a vector nelx by 1 of a unique random integer between 1 and p/2
% R = R.*2; %multiplies the random integers by 2 to create only even integers

R = randi(5,4,1);
Rstore(:,j) = R;
% R = [1,4,3,5]';
R = R'*2+((1:nelx)-1)*(nely+1)*2;


S(:,R(1:end)) = eye(nelx); %for all rows of S, a one is placed in the column specified by the random number of R


Fb = -1*ones(nelx,1);
Sp = zeros(p,nf); % Sp specifies the loading location
Sp(2*(1:nelx)*(nely+1),:) = eye(nelx); % put loads at the bottom of the beam

% Sptemp = Sp';
B = Sp(:,1)+Sp(:,2)+Sp(:,3)+Sp(:,4);%+Sp(:,5)+Sp(:,6)+Sp(:,7);%+Sp(:,8); %creates 1 vector with all Sp information
%NOTE: Must change B and C with every change in number of elements
B = B';
C = S(1,:)+S(2,:)+S(3,:)+S(4,:);%+S(5,:)+S(6,:)+S(7,:);%+S(8,:); %creates 1 vector with all S information
%NOTE: Must change B and C with every change in number of elements
strucS = zeros(p/nelx,nelx);
strucSp = zeros(p/nelx,nelx);
for i = 0:nelx-1
    strucS(:,i+1) = C(:,i*p/nelx+1:(i+1)*p/nelx); %creates physical structure resembling real set up of sensors
    strucSp(:,i+1) = B(:,i*p/nelx+1:(i+1)*p/nelx); %"..." of loads
    i = i+1;
end

[rowS,colS] = find(strucS); %finds location of sensors in the structure
[rowSp,colSp] = find(strucSp);%finds location of loads in the structure
rowS = rowS/2; %convert each node having 2 DOF to a single point
rowSp = rowSp/2;

Scoords = zeros(nelx,2);
Spcoords = zeros(nelx,2);

Scoords = [rowS colS]; %coordinates of the sensors
Spcoords = [rowSp colSp];%coordinates of the loads

Distance = pdist2(Scoords, Spcoords); %calculates distance from each sensor to each load
%ShortDist = sum(min(Distance)); %finds shortest sensor distance from each load
AvgShortDist(1,j) = mean(min(Distance')); %min goes by smallest in each row for this
% AvgShortDist = mean(ShortDist); %finds average of the shortest distances from each sensor to a load

% define material density
xPhys = rand(nely,nelx);
xPhys = xPhys/sum(xPhys(:))*volfrac*nelx*nely;
% xPhys = volfrac*ones(nely,nelx); 
% xPhys(2:end-1,2:end-1)=0;
% xPhys = [1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1];

% define the structure
% element-wise stiffness matrix for a quadrilateral element (square in shape)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]; %stiffness matrix
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

% element-wise consistent mass matrix for a quadrilateral element (square in shape)
ME = (4*eye(8) + [zeros(4),eye(4);eye(4),zeros(4)] + repmat(kron([0 1;1 0],2*eye(2)),2,2))/9;

% element-to-global assembly
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
sM = reshape(ME(:)*(Mmin+xPhys(:)'*(M0-Mmin)),64*nelx*nely,1);
M = sparse(iK,jK,sM); M = (M+M')/2;

Kb = K(freedofs,freedofs);
Mb = M(freedofs,freedofs);

D = S*inv(Mb)*Sp;
DReal(1,j) = trace(inv(D'*D+1e-6*eye(4))); %what does the eye(x) mean?
C1 = -S*inv(Mb)*Kb;
C = [zeros(num_observer,p), C1];
% M2 follows corollary 5 eq 35
% M1 follows eq 13b
j = j+1;

end
% LargeD = DReal>1.5e6; % Makes all values of DReal larger than 1.5e6=1, others=0
% Outliers = LargeD.*AvgShortDist; % combines LargeD and AvgShortDist
% Outliers(Outliers==0)=[]; % deletes entries equal to 0
% Outliers = sort(Outliers); % sorts data numerically
% sizeOutliers = numel(Outliers); % finds number of elements with DReal>1.5e6
% space = linspace(1,sizeOutliers, sizeOutliers); % creates linspace with sizeOutliers
% plot(space,Outliers)
% AvgOutlier = mean(Outliers); %finds average distance for a D value larger than 1.5e6


%DReal(DReal>1e6) = -1e6; %makes inf values equal to -1e6

plot(AvgShortDist,DReal,'.','MarkerSize',20)
xlabel('Average Distance Between Sensor and Input');
ylabel('Measure of Input Estimation Error');

