function [dT,varargout]=FLS(keq,Cs,r,H,D,q,varargin)
%
%% Description :
% The function computes the temperature at a point situated at any distance from one or many vertical finite line-sources.
% FLSM evaluates the response function of the model at selected times, interpolates at all evaluation times
% and makes a convolution in the frequency domain to compute the temperature.
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% FLSM Version 2.1.0 (March 21, 2013)
% Compatible with Matlab 7.14.0.739 (R2012a)
%
%% Reference :
% Marcotte, D. & Pasquier, P., 2014. Unit-response function for ground heat exchanger with parallel, series or mixed borehole arrangement. Renewable Energy, 68, pp.14–24.
%
%% References used to develop the code:
% 1) Claesson, J. & Javed, S., 2011. An analytical method to calculate borehole fluid temperatures for timescales from minutes to decades. In ASHRAE Annual Conference. Montréal, Canada, p. 10.
% 2) Marcotte, D. & Pasquier, P., 2008. Fast fluid and ground temperature computation for geothermal ground-loop heat exchanger systems. Geothermics, 37(6), pp.651–665.
%
%% Syntax :
% [dT,g_fft]=FLSM(keq,Cs,r,H,D,q,Method,g_fft)
%
%% Input and output variables :
%
% keq :        Equivalent ground thermal conductivity (W/mK)
% Cs  :        Volumetric heat capacity  (Cs=Cp*rho) (J/Km^3)
% r     :        Matrix of distances between the evaluation point and each of the nr sources [nr x nr] (m).
% H    :        Line-source length (m)
% D    :        Depth of the line-source (m)
% q     :        Unitary heating or cooling power of each line-source [nt x nr].
%                 The first column is the vector of time (s) where a step (of power) stops.  The steps are of equal duration.
%                 The following nr columns are the unitary heating (positive) or cooling (negative) power (W/m) of the nr lines-sources.
% dT  :         Temperature variation at a distance r from each line-source under the impulsion of signals q [nt x nr x nr] (K).
%                 The first dimension of dT corresponds to the time. The second dimension is associated to the heat sources and the third dimension corresponds to the evaluation distances.
%                 Temperature variation caused by nr line-sources is obtained by spatial superposition with the command squeeze(sum(dT,2));
% g_fft :       Optional - Spectrum of the response function [nt x nr x nr] (complex).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.0 -  Compute vector f and its spectrum
nt=size(q,1);
nr=size(r,2);
nr2=size(q,2)-1;
ti=q(:,1,1);

f=q(:,2:end)-[zeros(1,nr2);q(1:end-1,2:end)];
f_fft=fft(reshape(f,[nt 1 nr2]),2*nt,1); % Pad with 2*nt zeros

%% 2.0 - Select some evaluation time with a geometric progression
if nargin==6 % if g_fft is not provided construct function g
    if length(ti)>=10 % Subsample g if there is enough time steps
        id=2.^(0:log(length(ti))/log(2));
        if id(end)~=nt
            id=[id,nt];
        end
        id2=nan(length(id(1:end-1))*2,1);
        
        id2(1:2)=id(1:2);% Add the initial points
        id2(4:2:end)=id(3:end);
        id2(3:2:end-1)=diff(id(2:end))/2+id(2:end-1);
        id2=id2(mod(round(id2),id2)==0);
        
        t=ti(round(id2));
    else
        t=ti;
    end
elseif nargin==7 % Jump directly to convolution if g_fft is provided.
    g_fft=varargin{1};
end

%% 3.0 - Compute vector g and its spectrum
if nargin==6 % if g_fft is not provided construct function g
    r=reshape(r,[nr*nr2 1]); %Transform r
    [rr,~,Ind_c] = unique(r);    % Select a subset of r to reduce computation time
    
    parfor i=1:numel(t)     % Numerical quadrature of the formula published by Cleasson and Javed (2011)
        Temp=quadv(@(s)(exp(-rr.^2*s.^2)./(H*s.^2)).*(2*((H*s).*erf(H*s)-1/sqrt(pi)*(1-exp(-(H*s).^2)))+2*((H*s+2*D*s).*erf(H*s+2*D*s)-1/sqrt(pi)*(1-exp(-(H*s+2*D*s).^2)))-((2*H*s+2*D*s).*erf(2*H*s+2*D*s)-1/sqrt(pi)*(1-exp(-(2*H*s+2*D*s).^2)))-((2*D*s).*erf(2*D*s)-1/sqrt(pi)*(1-exp(-(2*D*s).^2)))),1./sqrt(4*keq/Cs.*t(i)),1e4);
        Integ(i,:,:)=reshape(Temp(Ind_c),[1 nr nr2]);
    end
    g=1/(4*pi*keq)*Integ;
    
    if nt~=numel(t) % Interpolate at each evaluation time by cubic spline
        gg=nan(nt,nr,nr2);
        for i=1:nr2
            gg(:,i,:)=spline(t,g(:,:,i)',ti)';
        end
    else
        gg=g;
    end
    g_fft=fft(gg,2*nt,1); % Pad with 2*nt zeros
    
    varargout{1}=g_fft;
    varargout{2}=nan;
end

%% 4.0 - Compute temperature variations in the spectral domain
dT=nan(2*nt,nr2,nr2);
clear f g gg q ti % free some memory space
for i=1:nr2
    dT(:,:,i)=ifft(f_fft.*g_fft(:,i,:),[],1);
end

dT=dT(1:nt,:,:);

end % End of function