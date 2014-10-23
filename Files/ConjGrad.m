function [x,f,g,nFoncObj,nLSearch]=ConjGrad(FoncObj,Options,varargin)
%
%% Description:
% ConjGrad determines the vector x minimizing the objective function F using the conjugate gradient method.
% Search along the descent direction is performed using the method of the golden section. \\
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% Version 1.5 (May 5,. 2014)
% Compatible with Matlab 8.2.0.701 (R2013b)
%.
%% Syntax :
% [x,f,nObjFun,nLSearch] =ConjGrad(@TRT_ObjFunc,Options,Para,Data,mean(Data.T_Exp(:,id),2),Inversion);
%
%% Input and output variables :
%
%  See  TRT_SInterp for variable description.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global mae_i f0_i
mae_i=inf;

if strcmp(varargin{1}.Options.Source,'FLSM')
    Xb=varargin{4}.FLSM.Xb(varargin{4}.FLSM.id,:);
    x=varargin{4}.FLSM.X0(varargin{4}.FLSM.id,evalin('caller','i'));
elseif strcmp(varargin{1}.Options.Source,'TRCM')
    Xb=varargin{4}.TRCM.Xb(varargin{4}.TRCM.id,:);
    x=varargin{4}.TRCM.X0(varargin{4}.TRCM.id,evalin('caller','i'));
end

[F(1),g] = FoncObj(x,varargin{:});
nFoncObj=numel(x)+1;
nLSearch=0;

optim=optimset('fminbnd');

for i = 1:Options.max_Iter
    if  mae_i<=Options.err_crit % If the first iteration converges
        f=f0_i;
        break;
    end
    
    % Evaluation of the descent direction
    if i==1   % Steepest descent for the first iteration
        d=-g;
    else      % Descent given by Hestenes-Stiefel for subsequent iterations
        d=-g+(g'*(g-g_o))/((g-g_o)'*d)*d;
        if sum(isnan(d))>0
            d=-g;
        end
    end
    g_o=g;
    
    % Limit the maximum value of t to keep the parameters within the boundaries specified by the user.
    tmax=([(Xb(:,1)-x)./d;(Xb(:,2)-x)./d]);
    id=tmax>=0 & not(isinf(tmax));
    tmax=min(tmax(id));
    if isempty(tmax)
        tmax=0.15;
    end
    
    % Line search
    optim=optimset(optim,'TolX',Options.TolRel_t*tmax,'OutputFcn', @stop_TRT_ObjFunc,'display','off');%,'PlotFcns',@optimplotfval
    [t,F(i+1),~,output]=fminbnd(FoncObj,0,tmax,optim,d,x,varargin{:});
    nFoncObj=nFoncObj+output.funcCount;
    nLSearch=nLSearch+1;
    
    x=x+t*d;
   
    % Check if stopping criteria is met or if a minimum is found.
    stop=(nFoncObj>=Options.max_Iter ||...
        mae_i<=Options.err_crit || ...
        (mae_i>=2*Options.err_crit && norm(d)<=0.4) || ...
        abs(F(i+1)-F(i))/abs(F(i))<Options.TolRel_F);
    
    if stop % Stopping criteria is met
        f=f0_i;
        break;
    else    % Compute F and G and loop again
        [F(i+1),g]=FoncObj(x,varargin{:});
        nFoncObj=nFoncObj+numel(x)+1;
    end
end
end