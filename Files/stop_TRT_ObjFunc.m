function stop = stop_TRT_ObjFunc(x,optimValues,state,varargin)
%
%% Description :
%
% stop_TRT_ObjFunc determines if the stopping criteria of an inversion is met and send a flag to stop the line search.
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

global mae_i
stop = false;
if strcmp(varargin{3}.Options.Source,'FLSM') && not(strcmp(state,'init'))
    stop=logical(mae_i<=varargin{6}.FLSM.SC_Inv);
elseif strcmp(varargin{3}.Options.Source,'TRCM') && not(strcmp(state,'init'))
    stop=logical(mae_i<=varargin{6}.TRCM.SC_Inv);
end
end