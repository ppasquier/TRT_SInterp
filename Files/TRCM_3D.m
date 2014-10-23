function [T,R,C,T0,t,Para]=TRCM_3D(Para)
%
%% Description :
% TRCM_3D generates a network of thermal resistances and capacites and integrates the resulting system of ODE.
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% TRCM_3D Version 2.5 (december 19, 2013)
% Compatible with Matlab 8.2.0.701 (R2013b)
%
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
% [T,R,C,T0,t,Para]=TRCM_3D(Para)
%
%% Input and output variables :
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.0 - Initialization of some variables

A=pi*(8*Para.r_b)^2;
Para.v=Para.Vdot_Fluid/(pi*Para.r_i^2);

nf=Para.n.nf;
np=Para.n.np;
ng=Para.n.ng;
ngg=Para.n.ngg;
nc=Para.n.nc;
ns=Para.n.ns;
nz=Para.n.nz;
nzs=Para.n.nzs;

n=2*(nf+np+ng)+ngg+nc+ns;

cg=Para.Cg;
cf=Para.Cp_Fluid*Para.rho_Fluid;
cp=Para.Cp_Pipe*Para.rho_Pipe;
cs=Para.Cs;

if Para.L==0
    Para.e=0.05;%[m]
    Para.k_Casing=Para.keq;
    Para.Cp_Casing=Para.Cp_Soil;
    Para.rho_Casing=Para.rho_Soil;
    
    cc=Para.Cs;
else
    cc=Para.Cp_Casing*Para.rho_Casing;
end

%% 2.0 - Compute the resistances of a layer
Para.Rs_z=(Para.dz/3)/(Para.keq*A);
Para.Cs_z=(Para.dz/3)*A*Para.rho_Soil*Para.Cp_Soil;

Rf=Para.Rf;
Rp=Para.Rp*Para.n_Pipe;
Ra=Para.R.MP_Ra;
Rb=Para.R.MP_Rb;
Rar=Ra-2*(Rf+Rp);
Rg=2*Rb-(Rf+Rp);
Rgg=2*Rg*((Rar))/(2*Rg-(Rar));

R=nan(2,n,2);
for j=1:2
    R(j,[1:nf,nf+np+ng+2:2*nf+np+ng+1],1)=Rf/nf;
    R(j,[nf+1:nf+np,2*nf+np+ng+2:2*nf+2*np+ng+1],1)=Rp/np;
    R(j,[nf+np+1:nf+np+ng,2*nf+2*np+ng+2:2*nf+2*np+2*ng+1],1)=Rg/ng;
    R(j,[nf+np+ng+1,2*nf+2*np+2*ng+2:2*nf+2*np+2*ng+ngg],1)=Rgg/ngg;
    if j==1
        Rs(1)=log(Para.r_Soil./(Para.r_b+Para.e))/(2*pi*Para.keq)/ns;
        Rc=log((Para.r_b+Para.e)./Para.r_b)/(2*pi*Para.k_Casing)/nc;
        R(j,2*nf+2*np+2*ng+ngg+1:2*nf+2*np+2*ng+ngg+nc,1)=Rc;
        R(j,2*nf+2*np+2*ng+ngg+nc+1:n,1)=Rs(1);
    else
        Rs(2)=log(Para.r_Soil./(Para.r_b))/(2*pi*Para.keq)/(ns+nc);
        R(j,2*nf+2*np+2*ng+ngg+1:n,1)=Rs(2);
    end
end

R=permute(R,[2 3 1]);
R=R(:,1,:);

%% 3.0  - Compute the volume of each node
rf=(1e-6);
V=nan(2,2*(nf+np+ng)+ngg+nc+ns);
kf_eq=log(Para.r_i/rf)/(2*pi*Rf);

vgg=2*Para.r_o*2*Para.D-(pi*Para.r_o^2);

rg_eq=sqrt((2*pi*Para.r_o^2+vgg)/pi);
kg_eq=log(Para.r_b/rg_eq)/(2*pi*Rg);% 2*Rb

% Volume of fluid
r=nan(1,nf);
vf=nan(1,nf);
r(1)=rf;
for i=1:nf
    r(i+1)=r(i)*(exp(2*pi*Rf/nf*kf_eq));
    vf(i)=pi*(r(i+1)^2-r(i)^2);
end
Para.r=[];
Para.r=[Para.r, r];

% Volume of pipe
r=nan(1,np);
vp=nan(1,np);
r(1)=Para.r_i;
for i=1:np
    r(i+1)=r(i)*(exp(2*pi*Rp/np*Para.k_Pipe));
    vp(i)=pi*(r(i+1)^2-r(i)^2);
end
Para.r=[Para.r, r];

% Volume of outer grout
r=nan(1,ng);
vg=nan(1,ng);
r(1)=rg_eq;
for i=1:ng
    r(i+1)=r(i)*(exp(2*pi*Rg/ng*kg_eq));
    vg(i)=pi*(r(i+1)^2-r(i)^2)/2;
end
Para.r=[Para.r, r];

% Volume of inner grout - Rectangle
vgg=vgg/ngg*ones(1,ngg);

for j=1:2
    V(j,[1:nf,nf+np+ng+2:2*nf+np+ng+1])=[vf,vf];
    V(j,[nf+1:nf+np,2*nf+np+ng+2:2*nf+2*np+ng+1])=[vp,vp];
    V(j,[nf+np+1:nf+np+ng,2*nf+2*np+ng+2:2*nf+2*np+2*ng+1])=[vg,vg];
    V(j,[nf+np+ng+1,2*nf+2*np+2*ng+2:2*nf+2*np+2*ng+ngg])=vgg;
    
    if j==1
        % Volume of casing
        r=nan(1,nc);
        vc=nan(1,nc);
        r(1)=Para.r_b;
        for i=1:nc
            r(i+1)=r(i)*(exp(2*pi*Rc*Para.k_Casing));
            vc(i)=pi*(r(i+1)^2-r(i)^2);
        end
        
        Para.r=[Para.r(1:end-1), r];
        
        % Volume of soil
        r=nan(1,ns);
        vs=nan(1,ns);
        r(1)=Para.r_b+Para.e;
        for i=1:ns
            if i==ns
                r(i+1)=Para.r_Soil;
            else
            r(i+1)=r(i)*(exp(2*pi*Rs(1)*Para.keq));
            end
            vs(i)=pi*(r(i+1)^2-r(i)^2);
        end
        
        Para.r=[Para.r, r(2:end)];
        
        V(j,2*nf+2*np+2*ng+ngg+1:2*(nf+np+ng)+ngg+nc)=vc;
        V(j,2*nf+2*np+2*ng+ngg+nc+1:2*(nf+np+ng)+ngg+nc+ns)=vs;
    else
        % Volume of soil
        r=nan(1,ns);
        vs=nan(1,ns);
        r(1)=Para.r_b;
        for i=1:nc+ns
            r(i+1)=r(i)*(exp(2*pi*Rs(2)*Para.keq));
            vs(i)=pi*(r(i+1)^2-r(i)^2);
        end
        V(j,2*nf+2*np+2*ng+ngg+1:2*(nf+np+ng)+ngg+nc+ns)=vs;
    end
end

%% 4.0 - Compute capacity of each node
C=nan(2,n,2);

for j=1:2
    C(j,[1,nf+np+ng+2],1)=cf*V(j,[1,nf+np+ng+2]);
    C(j,[nf+1,2*nf+np+ng+2],1)=cp*V(j,[nf+1,2*nf+np+ng+2])/2;
    C(j,[nf+2:nf+np,2*nf+np+ng+3:2*nf+2*np+ng+1],1)=cp*(V(j,[nf+1:nf+np-1,2*nf+np+ng+2:2*nf+2*np+ng])/2+V(j,[nf+2:nf+np,2*nf+np+ng+3:2*nf+2*np+ng+1])/2);
    C(j,[nf+np+1,2*nf+2*np+ng+2],1)=cp*V(j,[nf+np,2*nf+2*np+ng+1])/2+cg*V(j,[nf+np+1,2*nf+2*np+ng+2])/2+cg*V(j,[nf+np+ng+1,2*nf+2*np+2*ng+ngg])/2;
    C(j,[nf+np+2:nf+np+ng,2*nf+2*np+ng+3:2*nf+2*np+2*ng+1],1)=cg*(V(j,[nf+np+1:nf+np+ng-1,2*nf+2*np+ng+2:2*nf+2*np+2*ng])/2+V(j,[nf+np+2:nf+np+ng,2*nf+2*np+ng+3:2*nf+2*np+2*ng+1])/2);
    C(j,[2*nf+2*np+2*ng+2],1)=cg*(V(j,[nf+np+ng+1])/2+V(j,[2*nf+2*np+2*ng+2])/2);
    C(j,[2*nf+2*np+2*ng+3:2*nf+2*np+2*ng+ngg],1)=cg*(V(j,[2*nf+2*np+2*ng+2:2*nf+2*np+2*ng+ngg-1])/2+V(j,[2*nf+2*np+2*ng+3:2*nf+2*np+2*ng+ngg])/2);
    C(j,[nf+np+ng+1],1)=cg*(V(j,[nf+np+ng])/2+V(j,[2*nf+2*np+2*ng+1])/2)+cc*V(j,[2*nf+2*np+2*ng+ngg+1])/2;
    C(j,[2*nf+2*np+2*ng+ngg+1:2*(nf+np+ng)+ngg+nc-1],1)=cc*(V(j,[2*nf+2*np+2*ng+ngg+1:2*(nf+np+ng)+ngg+nc-1])/2+V(j,[2*nf+2*np+2*ng+ngg+2:2*(nf+np+ng)+ngg+nc])/2);
    C(j,[2*(nf+np+ng)+ngg+nc],1)=cc*(V(j,[2*(nf+np+ng)+ngg+nc])/2)+cs*V(j,[2*nf+2*np+2*ng+ngg+nc+1])/2;
    C(j,[2*nf+2*np+2*ng+ngg+nc+1:2*(nf+np+ng)+ngg+nc+ns-1],1)=cs*(V(j,[2*nf+2*np+2*ng+ngg+nc+1:2*(nf+np+ng)+ngg+nc+ns-1])/2+V(j,[2*nf+2*np+2*ng+ngg+nc+2:2*(nf+np+ng)+ngg+nc+ns])/2);
    C(j,[2*(nf+np+ng)+ngg+nc+ns],1)=cs*(V(j,[2*(nf+np+ng)+ngg+nc+ns])/2);
end

C=permute(C,[2 3 1]);
C=C(:,1,:);

%% 5.0 - Prepare initial conditions
T0=Para.Tg*ones(1,n*(nz)+nzs);
T0(1)=Para.Tg+Para.dT;

%% 6.0 - Integration of ODE system
options=odeset('RelTol',1e-3,'abstol',1e-6,'stats','off','vectorize','on');
[t,T]=ode15s(@Tdot_TRCM_3D,[Para.t]',T0,options,R,C,Para);

Para.C=C;
Para.R.R_TRCM=R;
Para=orderfields(Para);
