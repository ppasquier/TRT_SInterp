function [T_TRCM,g_TRCM,g_Data,Para]=TRCM_Spectral(Para,Data)

%% Description :
%
% TRCM_Spectral convolves, in the spectral domain, an impulse function f and a transfer function g computed by a quasi-3D thermal resistance and capacity model.
% TRCM.
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
%   [T_TRCM,g_TRCM,g_Data,Para]=TRCM_Spectral(Para,Data)
%
%% Input variables :
%
% Para: Structure variable containing the parameters of the BHE
%            r_b : Borehole radius (m)
%            r_i :  Inner pipe radius (m)
%            r_o : Outer pipe radius (m)
%            r_Soil : Outer radius of the model (m)
%            D:     Distance between the center of the borehole and the center of a pipe (m)
%            H:     Borehole length (m)
%           Cp_Fluid : Specific heat of the fluid (J/kg K)
%           Cp_Grout : Specific heat of the filling material (J/kg K)
%           Cp_Pipe : Specific heat of the pipe (J/kg K)
%           Cp_Soil : Specific heat of the soil (J/kg K)
%           rho_Fluid : Density of the fluid (kg/m^3)
%           rho_Grout : Density of the filling material (kg/m^3)
%           rho_Pipe : Density of the pipe (kg/m^3)
%           rho_Soil : Density of the soil (kg/m^3)
%           k_Grout : Thermal conductivity of the filling material (W/mK)
%           k_Pipe : Thermal conductivity of the pipe (W/mK)
%           k_Soil : Thermal conductivity of the soil (W/mK)
%           Rf : Convective fluid resistance (mK/W)
%           Tg:  Initial ground temperature (o^C)
%           Vdot_Fluid : Circulation flow rate (m^3/s)
%
%   See  TRT_SInterp for a description of the variable Data
%
%% 1.0 - Variable Initialization

Para.n_Pipe=2;
Para.Cg=Para.rho_Grout*Para.Cp_Grout;
Para.Cs=Para.rho_Soil*Para.Cp_Soil;
Para.keq=Para.k_Soil;

Para.n.nf=1;
Para.n.np=6;
Para.n.ng=6;
Para.n.ngg=6; % Intersection supérieure - Borehole wall
Para.n.nc=6;
Para.n.ns=4*6;
Para.n.nz=max(round(Para.H/4),50);
Para.n.nzs=10;
n=Para.n;
nn=2*(n.nf+n.np+n.ng)+n.ngg+n.nc+n.ns;

Para.dz=Para.H/(n.nz-1);
Para.z=[0:Para.dz:Para.H];

Para.v_Fluid=Para.Vdot_Fluid/(pi*Para.r_i.^2);
if ~isfield(Para,'Rf')
    Para.Rf=1e-9;
end

%% 2.0 - Compute equivalent borehole resistance

Para.Rp=log(Para.r_o/Para.r_i)/(2*pi*Para.k_Pipe*Para.n_Pipe);
[Para.R]=BoreholeR(Para.D,Para.k_Grout,Para.k_Pipe,Para.k_Soil,Para.r_b,Para.r_i,Para.r_o,Para.n_Pipe,Para.H,nan,Para.Cp_Fluid,Para.Rf,Para.Rp);

%% 3.0 - Compute minimum time step

Para.dt=min([Data.t(2)-Data.t(1) round((2*Para.H/Para.v_Fluid)/12)]);
Para.t=[Para.dt:Para.dt:Data.t(end),Data.t(end)]';
Para.t=unique(Para.t);
Para.nt=numel(Para.t);

Para.delta_T=interp1(Data.t,Data.delta_T,Para.t,'linear',Data.delta_T(1));

%% 4.0 - Compute transfer function g

Para.dT=2;
[g_TRCM,~,~,~,~,Para]=TRCM_3D(Para);
g_TRCM=(g_TRCM-Para.Tg)/Para.dT;

%% 5.0 - Interpolate transfer function g

g=nan(Para.nt,size(Data.T_Probe_Location,2));
[ind]=unique(Data.T_Probe_Location(2,:));

for i=1:size(ind,2)
    if ind(i)==1
        T_z=g_TRCM(:,2*(n.nf+n.np+n.ng+1)+0:nn:(n.nz)*nn);
    elseif ind(i)==2
        T_z=g_TRCM(:,(1)+0:nn:(n.nz)*nn);
    elseif ind(i)==3
        T_z=g_TRCM(:,(n.nf+n.np+n.ng+1)+0:nn:(n.nz)*nn);
    elseif ind(i)==4
        T_z=g_TRCM(:,(2*n.nf+2*n.np+2*n.ng+n.ngg)+0:nn:(n.nz)*nn);
    elseif ind(i)==5
        T_z=g_TRCM(:,(n.nf+n.np+n.ng+2)+0:nn:(n.nz)*nn);
    elseif ind(i)==6
        T_z=g_TRCM(:,(2*n.nf+2*n.np+2*n.ng+n.ngg)+0:nn:(n.nz)*nn);
    end
    
    F=griddedInterpolant({Para.t,Para.z},T_z);
    
    ind2=find(ind(i)==Data.T_Probe_Location(2,:));
    for j=1:size(ind2,2)
        g(:,ind2(j))=F(Para.t,Data.T_Probe_Location(1,ind2(j))*ones(Para.nt,1));
    end
end

g_Data=g;

% Tz_TRCM=[T(:,1:nn:end)';fliplr(T(:,(n.nf+n.np+n.ng+2):nn:end))'];
 
%% 0.0 - Computation of the input function f
f=[Para.delta_T(1);diff(Para.delta_T)];

%% 0.0 - Convolution in the spectral domain of the zero padded functions f and g 
T_Num=ifft(bsxfun(@times,fft(f,2*numel(f)),fft(g,2*size(g,1))));
T_Num=T_Num(1:numel(f),:)+Para.Tg;

%% 0.0 - Interpolation at the time of the data
T_TRCM=interp1(Para.t,T_Num,Data.t,'linear','extrap');
end

function [R]=BoreholeResistance(D,kg,kp,keq,rb,ri,ro,npipe,H,m,Cp,Rf,Rp)
Rpn=npipe*[Rf+Rp];

dtheta=0.5;
theta=[dtheta:dtheta:360];
dr=1e-6;
dtheta=0.5;

if npipe==2
    xym=[D 0;-D 0];
    rp=[ro ro];
    Rpn=[Rpn Rpn];
    To=[35 20];
elseif npipe==4
    xym=[D 0;-D 0;0 D;0 -D];
    rp=[ro ro ro ro];
    Rpn=[Rpn Rpn Rpn Rpn];
    To=[35 20 35 20];
end

xg=(rb+dr)*cosd(theta)';
yg=(rb+dr)*sind(theta)';
[Taa,qaa]=Multipole_J(keq,kg,rb,rp,To,xym,[xg yg],Rpn,10);

dr=-dr;
xg=(rb+dr)*cosd(theta)';
yg=(rb+dr)*sind(theta)';
[Tbb,qbb]=Multipole_J(keq,kg,rb,rp,To,xym,[xg yg],Rpn,10);

Tb=mean(mean([Taa,Tbb],2));
Tb1=mean(Taa([1:90/dtheta (270+dtheta)/dtheta:360/dtheta]));
Tb2=mean(Taa([(90+dtheta)/dtheta:270/dtheta]));

R.MP_Rb_eff=(mean(To)-Tb)/sum(mean([qaa;qbb],1));

if npipe==2
    R1=(To(1)+To(2)-(Tb1+Tb2))/(qaa(1)+qaa(2));
    R12=2*(To(1)-To(2))/(qaa(1)-qaa(2)-(To(1)-To(2)-(Tb1-Tb2))/R1);
    R12b=2*(Tb1-Tb2)/(qaa(1)-qaa(2)-(To(1)-To(2)-(Tb1-Tb2))/R1);
    
    Rb=R1/2;
    Ra=4*R12*Rb/(4*Rb+R12);
    
    R.MP_Ra=Ra;
    R.MP_Rb=Rb;
    R.MP_R12=R12;
end
end
