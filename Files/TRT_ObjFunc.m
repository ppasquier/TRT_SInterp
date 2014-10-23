function [varargout]=TRT_ObjFunc(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description :
%
% TRT_ObjFunc computes the objective function F and its gradient G by finite differences for the finite line-source model and thermal resistance and capacity model.
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% Version 2.8 (May 5,. 2014)
% Compatible with Matlab 8.2.0.701 (R2013b)
%.
%% Input and output variables :
%
%  See  TRT_SInterp for variable description.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.0 - Variable Initialization

global mae_i f0_i

if nargin==5    % Evaluation of F and G
    [x,Para,Data,T_Experimental,Inversion]=deal(varargin{:});
else                % Evaluation of F only
    [t,d,x,Para,Data,T_Experimental,Inversion]=deal(varargin{:});
    x = x + t*d;
end

ind_inv=Inversion.cpt;

if strcmp(Para.Options.Source,'FLSM')
    ind=find(Inversion.FLSM.id);
else
    ind=find(Inversion.TRCM.id);
end

pas=0.0001;
dx=ones(6,6);
G=zeros(numel(x),1);
for i=1:numel(x)
    dx(ind(i),ind(i))=dx(ind(i),ind(i))+pas;
end

%% 2.0 - Evaluation of F with the FLS model
if strcmp(Para.Options.Source,'FLSM')   
    Rb=@(Para,ks,kg,D)(1/(4*pi*kg)*(log(Para.r_b/Para.r_o)+log(Para.r_b/2./D)+(kg-ks)/(kg+ks)*log((Para.r_b./D).^4./((Para.r_b./D).^4-1)))+Para.Rp);

    q=Data.Qdot/Para.H;
    D=mean(Inversion.Xb(5,:),2);
    X=Inversion.FLSM.X0(:,ind_inv);
    X(Inversion.FLSM.id)=x;       
        
    % Evaluation of the objective function F
    f0=X(4)+q.*Rb(Para,X(1),X(3),D)+FLS(X(1),X(2),Para.r_b,Para.H,0,[Data.t q]);
    mae=mean(mean(abs(T_Experimental(Para.id_F)-f0(Para.id_F))));
    F0=sqrt(mean(mean((T_Experimental(Para.id_F)-f0(Para.id_F)).^2)));
                    
    % Evaluation of G, the gradient of F
    if nargin<=5
        for i=1:sum(Inversion.FLSM.id)
            f=X(4)*dx(4,ind(i))+q.*Rb(Para,X(1)*dx(1,ind(i)),X(3)*dx(3,ind(i)),D)+FLS(X(1)*dx(1,ind(i)),X(2)*dx(2,ind(i)),Para.r_b,Para.H,0,[Data.t q]);
            F=sqrt(mean(mean((T_Experimental(Para.id_F)-f(Para.id_F)).^2)));
                   
            G(i)=-(F0-F)./(X(ind(i))*pas);
        end
    end
    
    %% 3.0 - Evaluation of F with TRCM
elseif strcmp(Para.Options.Source,'TRCM')
    
    X=Inversion.TRCM.X0(:,ind_inv);
    X(Inversion.TRCM.id)=x;
        
    Para.Vdot_Fluid= Data.Vdot;
    Para.k_Soil=X(1);
    Para.Cp_Soil=X(2)/1000;
    Para.k_Grout=X(3);
    Para.Cp_Grout=X(4)/1000;
    Para.D=X(5);
    Para.Tg=X(6);
    
    % Evaluation of the objective function F
    [f0,~,~,Para]=TRCM_Spectral(Para,Data);
    mae=mean(mean(abs(Data.T_Exp(Para.id_F,:)-f0(Para.id_F,:))));
    F0=sqrt(mean(mean((Data.T_Exp(Para.id_F,:)-f0(Para.id_F,:)).^2)));
    
    % Evaluation of G, the gradient of F
    if nargin<=5
        for i=1:sum(Inversion.TRCM.id)
            Para.k_Soil=X(1)*dx(1,ind(i));
            Para.Cp_Soil=X(2)*dx(2,ind(i))/1000;
            Para.k_Grout=X(3)*dx(3,ind(i));
            Para.Cp_Grout=X(4)*dx(4,ind(i))/1000;
            Para.D=X(5)*dx(5,ind(i));
            Para.Tg=X(6)*dx(6,ind(i));
            
            [f,~,~,~]=TRCM_Spectral(Para,Data);
            F=sqrt(mean(mean((Data.T_Exp(Para.id_F,:)-f(Para.id_F,:)).^2)));
            
            G(i)=-(F0-F)./(X(ind(i))*pas);
        end
    end
end

mae_i=mae;
f0_i=f0;

varargout{1}=F0;

if nargin<=5
    varargout{2}=G;
end
end