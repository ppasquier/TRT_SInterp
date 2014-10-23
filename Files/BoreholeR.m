function [R]=BoreholeR(D,kg,kp,keq,rb,ri,ro,npipe,H,m,Cp,Rf,Rp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description of the function :
% BoreholeR computes the equivalent borehole resistance with the multipole method of Bennet et al. (1987).
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% BoreholeR Version 1.1.0 (03-Mar-2011 10:03)
% Compatible avec Matlab 7.11.0.584 (R2010b)
%
%% References :
% 1) Hellström, G., 1991. Ground Heat Storage. Thermal Analysis of Duct Storage Systems. Part I Theory. University of Lund,  Sweden.
% 2) Bennet, J., Claesson, J. & Hellström, G., 1987. Multipole Method to Compute the Conductive Heat Flows to and between Pipes in a Composite
%     Cylinder, Note on Heat Transfer 3, Lund, Sweden: University of Lund, Department of Building Technology and Mathematical Physics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.0 - Use multipole method of Bennet et al.

Rpn=npipe*[Rf+Rp];

% Compute temperature for 720 points along the borehole wall
dtheta=0.5;
theta=[dtheta:dtheta:360];
dr=1e-6;
dtheta=0.5;

if npipe==2 % For 2 pipes
    xym=[D 0;-D 0];
    rp=[ro ro];
    Rpn=[Rpn Rpn];
    To=[35 20];
elseif npipe==4 % For 4 pipes
    xym=[D 0;-D 0;0 D;0 -D];
    rp=[ro ro ro ro];
    Rpn=[Rpn Rpn Rpn Rpn];
    To=[35 20 35 20];
end

% Compute at r+dr
xg=(rb+dr)*cosd(theta)';
yg=(rb+dr)*sind(theta)';
[Taa,qaa]=Multipole_J(keq,kg,rb,rp,To,xym,[xg yg],Rpn,10);

% Compute at r-dr
dr=-dr;
xg=(rb+dr)*cosd(theta)';
yg=(rb+dr)*sind(theta)';
[Tbb,qbb]=Multipole_J(keq,kg,rb,rp,To,xym,[xg yg],Rpn,10);

% Compute mean value
Tb=mean(mean([Taa,Tbb],2));
Tb1=mean(Taa([1:90/dtheta (270+dtheta)/dtheta:360/dtheta]));
Tb2=mean(Taa([(90+dtheta)/dtheta:270/dtheta]));

% Compute Rb
R.MP_Rb_eff=(mean(To)-Tb)/sum(mean([qaa;qbb],1));

% Compute internal resistances
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