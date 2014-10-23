function Tdot =Tdot_TRCM_3D(t,T,R,C,Para)

%% Description   :
% Tdot_TRCM_3D is an incredibly clean chunk of code full of very usefull comments and explanations used by TRCM_3D to compute a matrix of dT/dt for each node of the TRCM network. 
% 
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% Tdot_TRCM_3D Version 1.0.0 (December 19, 2012)
% Compatible with Matlab 7.14.0.739 (R2012a)
%.
%% Reference :
% Please, cite this work as :
% Pasquier, P. & Marcotte, D., 2013. Joint Use of Quasi-3D Response Model and Spectral Method to Simulate Borehole Heat Exchanger. Geothermics.
%
%% Additional information can be found in the following references :
% 1) Pasquier, P. & Marcotte, D., 2012. Short-term simulation of ground heat exchanger with an improved TRCM. Renewable Energy, 46, pp.92–99.
% 2) Bauer, D. et al., 2010. Thermal resistance and capacity models for borehole heat exchangers. International Journal of Energy Research, 35(4), pp.312–320.
% 3) Marcotte, D. & Pasquier, P., 2008. Fast fluid and ground temperature computation for geothermal ground-loop heat exchanger systems. Geothermics, 37(6), pp.651–665.
% 4) Bennet, J., Claesson, J. & Hellström, G., 1987. Multipole Method to Compute  the Conductive Heat Flows to and between Pipes in a Composite Cylinder, Lund, Sweden: University of Lund, Department of Building Technology and Mathematical Physics.
%
%% Syntax :
%Tdot =Tdot_TRCM_3D(t,T,R,C,Para)
%
%% Input and output variables :
% See TRCM 3D for a description of inputs and outputs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Initialization of some variables
% 1.1 - Number of nodes 
nf=Para.n.nf; 
np=Para.n.np;
ng=Para.n.ng;
ngg=Para.n.ngg;
nc=Para.n.nc;
ns=Para.n.ns;
nz=Para.n.nz;
nzs=Para.n.nzs;
n=2*(nf+np+ng)+ngg+nc+ns;

% 1.2 - Initialization of resistances
dT=Para.dT;
v=Para.v;
z=Para.z;
dz=Para.dz;
R_Loop=1e-6;
Rs_z=Para.Rs_z;
Cs_z=Para.Cs_z;
R_vertical=10;

R0=repmat(R/dz,[1 size(T,2) 1]);
C=repmat(C*dz,[1 size(T,2) 1]);

% 1.3 - Initialization of Tdot

Tdot=nan(size(T));

%% 2.0 - Construction of Tdot.
for k=0:nz-1 % Loop on the number of layers
    if k==0 && nz~=1
        R=0.5*R0;
    else
        R=R0;
    end
    
    if z(k+1)<=Para.L && Para.L~=0
        j=1;
    else
        j=2;
    end
    
    % 2.1 - Nodes of the downward pipe
    ind=1+k*n;
    if k==0 && nz~=1 % Surface
        Tdot(ind,:)=((T(ind+1,:)-T(ind,:))./R(ind,:,j)+v*0.5*C(ind-k*n,:,j).*(T(ind+n,:)-T(ind,:))/dz+(T(nf+np+ng+2,:)+dT-T(ind,:))./R_Loop)./(C(ind-k*n,:,j)); % Température variable
    elseif k~=0 && nz~=1 && nz~=2 && k~=(nz-1) % Intermediate layers
        Tdot(ind,:)=((T(ind+1,:)-T(ind,:))./R(nf,:,j)+v*C(ind-k*n,:,j).*(T(ind-n,:)-T(ind,:))/dz)./C(ind-k*n,:,j);
    elseif k==(nz-1) && nz~=1 % U-Loop.
        Tdot(ind,:)=((T(ind+1,:)-T(ind,:))./R(nf,:,j)+v*C(ind-k*n,:,j).*(T(ind-n,:)-T(ind,:))/dz)./(C(ind-k*n,:,j));
    end
    
    % 2.2 - Nodes of the U-Loop.
    if k==(nz-1) && nz~=1
        ind=nz*n+1;
        Tdot(ind,:)=((T(nf+np+ng+1+k*n,:)-T(ind,:))./Rs_z+(T(ind+1,:)-T(ind,:))./Rs_z)./Cs_z;
    end
    
    % 2.3 - Nodes of the BHE
    ind=[2:nf+np]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2:nf+np]-1,:,j)+(T(ind+1,:)-T(ind,:))./R([2:nf+np],:,j))./C(ind-k*n,:,j);
    
    ind=[nf+np+1]+k*n; % Left intersection
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R(nf+np,:,j)+(T(ind+1,:)-T(ind,:))./R(nf+np+1,:,j)+(T(2*(nf+np+ng+1)+k*n,:)-T(ind,:))./R(nf+np+ng+1,:,j))./C(ind-k*n,:,j);
    ind=[nf+np+2:nf+np+ng]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([nf+np+2:nf+np+ng]-1,:,j)+(T(ind+1,:)-T(ind,:))./R([nf+np+2:nf+np+ng],:,j))./C(ind-k*n,:,j);
    
    % 2.4 - Intersection - Borehole wall
    ind=[nf+np+ng+1]+k*n;
    if k==0 && strcmp(Para.Options.BC,'Adiabatic') && nz~=1 %Surface - Adiabatic
        Tdot(ind,:)=((T(nz*n+1,:)-T(ind,:))./Rs_z+(T(ind-1,:)-T(ind,:))./R([nf+np+ng+1]-1,:,j)+(T(2*nf+2*np+2*ng+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+1,:,j)+(T(2*nf+2*np+2*ng+ngg+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+ngg+1,:,j))./(Cs_z/2+C(ind-k*n,:,j));
    elseif k==0 && strcmp(Para.Options.BC,'TConstant') && nz~=1 %Surface - Constant temperature
        Tdot(ind,:)=((Para.Tg-T(ind,:))/R_vertical+(T(nz*n+1,:)-T(ind,:))./Rs_z+(T(ind-1,:)-T(ind,:))./R([nf+np+ng+1]-1,:,j)+(T(2*nf+2*np+2*ng+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+1,:,j)+(T(2*nf+2*np+2*ng+ngg+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+ngg+1,:,j))./(Cs_z/2+C(ind-k*n,:,j));
    elseif k==(nz-1) && nz~=1 % Last layer and link below BHE
        Tdot(ind,:)=((T(nz*n+1,:)-T(ind,:))./Rs_z+(T(ind-1,:)-T(ind,:))./R([nf+np+ng+1]-1,:,j)+(T(2*nf+2*np+2*ng+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+1,:,j)+(T(2*nf+2*np+2*ng+ngg+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+ngg+1,:,j))./(Cs_z/2+C(ind-k*n,:,j));
    else
        Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([nf+np+ng+1]-1,:,j)+(T(2*nf+2*np+2*ng+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+1,:,j)+(T(2*nf+2*np+2*ng+ngg+1+k*n,:)-T(ind,:))./R(2*nf+2*np+2*ng+ngg+1,:,j))./C(ind-k*n,:,j);
    end
    
    % 2.5 - Nodes of the upward pipe
    ind=nf+np+ng+2+k*n;
    if k==0 && nz~=1 %Surface
        Tdot(ind,:)=((T(ind+1,:)-T(ind,:))./R(nf+np+ng+2,:,j)+v*0.5*C(ind-k*n,:,j).*(T(ind+n,:)-T(ind,:))/dz)./(C(ind-k*n,:,j));
    elseif k~=0 && nz~=1 && nz~=2 && k~=(nz-1) % Intermediate layers
        Tdot(ind,:)=((T(ind+1,:)-T(ind,:))./R(nf+np+ng+2,:,j)+v*C(ind-k*n,:,j).*(T(ind+n,:)-T(ind,:))/dz)./C(ind-k*n,:,j);
    elseif k==(nz-1) && nz~=1 %U-Loop.
        Tdot(ind,:)= Tdot(1+k*n,:);
    end
    
    ind=[nf+np+ng+3:2*nf+2*np+ng+1]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([nf+np+ng+3:2*nf+2*np+ng+1]-1,:,j)+(T(ind+1,:)-T(ind,:))./R([nf+np+ng+3:2*nf+2*np+ng+1],:,j))./C(ind-k*n,:,j);
    
    % 2.6 - Right intersection
    ind=[2*nf+2*np+ng+2]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2*nf+2*np+ng+1],:,j)+(T(ind+1,:)-T(ind,:))./R([2*nf+2*np+ng+2],:,j)+(T(2*nf+2*np+2*ng+ngg+k*n,:)-T(ind,:))./R([2*nf+2*np+2*ng+ngg],:,j))./C(ind-k*n,:,j);
    ind=[2*nf+2*np+ng+3:2*(nf+np+ng)]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2*nf+2*np+ng+3:2*(nf+np+ng)]-1,:,j)+(T(ind+1,:)-T(ind,:))./R([2*nf+2*np+ng+3:2*(nf+np+ng)],:,j))./C(ind-k*n,:,j);
    ind=[2*(nf+np+ng)+1]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2*(nf+np+ng)],:,j)+(T((nf+np+ng)+1+k*n,:)-T(ind,:))./R([2*(nf+np+ng)+1],:,j))./C(ind-k*n,:,j);
    
    ind=[2*(nf+np+ng+1)]+k*n;
    Tdot(ind,:)=((T(nf+np+1+k*n,:)-T(ind,:))./R((nf+np+ng+1),:,j)+(T(ind+1,:)-T(ind,:))./R(2*(nf+np+ng+1),:,j))./C(ind-k*n,:,j);
    
    ind=[2*(nf+np+ng+1)+1:2*(nf+np+ng)+ngg-1]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2*(nf+np+ng+1)+1:2*(nf+np+ng)+ngg-1]-1,:,j)+(T(ind+1,:)-T(ind,:))./R([2*(nf+np+ng+1)+1:2*(nf+np+ng)+ngg-1],:,j))./C(ind-k*n,:,j);
    
    ind=[2*(nf+np+ng)+ngg]+k*n;
    Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2*(nf+np+ng)+ngg]-1,:,j)+(T((2*nf+2*np+ng)+2+k*n,:)-T(ind,:))./R([2*(nf+np+ng)+ngg],:,j))./C(ind-k*n,:,j);
    
    % 2.7 - Node of the casing and ground
    if (k==0 || k==(nz-1)) && strcmp(Para.Options.BC,'TConstant') && nz~=1 
        % First node of the casing
        ind=[2*(nf+np+ng)+ngg+1]+k*n;
        Tdot(ind,:)=((Para.Tg-T(ind,:))/R_vertical+(T(nf+np+ng+1+k*n,:)-T(ind,:))./R(2*(nf+np+ng)+ngg+1,:,j)+(T(ind+1,:)-T(ind,:))./R(2*(nf+np+ng)+ngg+2,:,j))./C(ind-k*n,:,j);
        
        % Casing + Ground
        ind=[2*(nf+np+ng+1)+ngg:2*(nf+np+ng)+ngg+nc+ns-1]+k*n;
        Tdot(ind,:)=((Para.Tg-T(ind,:))/R_vertical+(T(ind-1,:)-T(ind,:))./R([2*(nf+np+ng+1)+ngg:2*(nf+np+ng)+ngg+nc+ns-1],:,j)+(T(ind+1,:)-T(ind,:))./R([2*(nf+np+ng+1)+ngg+1:2*(nf+np+ng)+ngg+nc+ns],:,j))./C(ind-k*n,:,j);
        ind=(2*(nf+np+ng)+ngg+nc+ns)+k*n;
        Tdot(ind,:)=0;
        
        % Under the borehole
        ind=(nz*n+2):(nz*n+nzs-1);
        Tdot(ind,:)=0;
        
        % Last node of the ground
        ind=nz*n+nzs;
        Tdot(ind,:)=0;
    else
        % First node of the casing
        ind=[2*(nf+np+ng)+ngg+1]+k*n;
        Tdot(ind,:)=((T(nf+np+ng+1+k*n,:)-T(ind,:))./R(2*(nf+np+ng)+ngg+1,:,j)+(T(ind+1,:)-T(ind,:))./R(2*(nf+np+ng)+ngg+2,:,j))./C(ind-k*n,:,j);
        
        % Casing + Ground
        ind=[2*(nf+np+ng+1)+ngg:2*(nf+np+ng)+ngg+nc+ns-1]+k*n;
        Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./R([2*(nf+np+ng+1)+ngg:2*(nf+np+ng)+ngg+nc+ns-1],:,j)+(T(ind+1,:)-T(ind,:))./R([2*(nf+np+ng+1)+ngg+1:2*(nf+np+ng)+ngg+nc+ns],:,j))./C(ind-k*n,:,j);
        ind=(2*(nf+np+ng)+ngg+nc+ns)+k*n;
        Tdot(ind,:)=0;
        
        % Under the borehole
        ind=(nz*n+2):(nz*n+nzs-1);
        Tdot(ind,:)=((T(ind-1,:)-T(ind,:))./Rs_z+(T(ind+1,:)-T(ind,:))./Rs_z)./Cs_z;
        
        % Last node of the ground
        ind=nz*n+nzs;
        Tdot(ind,:)=0;
    end
end