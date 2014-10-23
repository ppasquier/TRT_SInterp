function [X,T,Stats]=TRT_SInterp(Data,Para,Inversion)
%
%% Description of the function :
% TRT_SInterp finds the combination of parameters minimizing the calibration error between the temperatures measured during a thermal response test
% and the temperatures predicted by a finite line-source model (FLSM) and/or a thermal resistance and capacity model (TRCM).  When using the TRCM,
% the program allows direct integration of temperature measurements made at different depths in the heat carrier fluid and borehole heat exchanger.
% The FLSM however uses only the mean fluid temperature measured by the probes located in the circulation pipes to evaluate the objective function of
% the optimization problem. To prevent overfitting of the models, the optimization process is stopped  as soon as the mean absolute calibration error falls
% below a stopping criteria specified by the user.  This criteria should be comparable to the uncertainty of the measurement probes.
% The program allows realization of  stochastic inversions to address the non-uniqueness of the inverse problem. It is however possible to do deterministic
% inversions by limiting the maximum number of realization of the Monte-Carlo experiment to one.
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% Version 1.1.0 (10-July-2014)
% Compatible with Matlab 7.14.0.739 (R2013b)
%
%% Reference :
% Please, cite this work as
% Pasquier, P., 2014. Stochastic interpretation of thermal response test with TRT-SInterp. Submitted to Computer & Geosciences.
%
%% Additional information can be found in the following references :
% 1) Pasquier, P. & Marcotte, D., 2013. Joint Use of Quasi-3D Response Model and Spectral Method to Simulate Borehole Heat Exchanger. Geothermics.
% 2) Pasquier, P. & Marcotte, D., 2012. Short-term simulation of ground heat exchanger with an improved TRCM. Renewable Energy, 46, pp.92–99.
% 3) Claesson, J. & Javed, S., 2011. An analytical method to calculate borehole fluid temperatures for timescales from minutes to decades. In ASHRAE Annual Conference. Montréal, Canada, p. 10.
% 4) Marcotte, D. & Pasquier, P. 2008.  On the estimation of thermal resistance in borehole conductivity test. Renewable Energy, vol. 33, p. 2407-2415. doi:10.1016/j.renene.2008.01.021
% 5) Marcotte, D. & Pasquier, P., 2008. Fast fluid and ground temperature computation for geothermal ground-loop heat exchanger systems. Geothermics, 37(6), pp.651–665.
% 6) Hellström, G., 1991. Ground Heat Storage. Thermal Analysis of Duct Storage Systems. Part I Theory. University of Lund,  Sweden.
% 7) Bennet, J., Claesson, J. & Hellström, G., 1987. Multipole Method to Compute  the Conductive Heat Flows to and between Pipes in a Composite Cylinder, Lund, Sweden: University of Lund, Department of Building Technology and Mathematical Physics.
%
%% Syntax :
%  [X,T,Stats]=TRT_SInterp(Data,Para,Inversion)
%
%% Input variables:
% Data: Structure array containing the following fields
%       dt:                             Scalar containing the time interval (s) between the measurements.  
%       T_Probe_Location:  2 x n_p matrix containing the location of the probes in the BHE. The first row is the installation depth (m) of each probe.  
%                                       The second row of the matrix is an integer flag ranging from 1 to 6  and specifying the horizontal location of the probe in the BHE. 
%                                                                              Tin   Tout
%                                                                            /¯¯¯¯¯¯¯¯\
%                                                                           (   (x)     (.)  )
%                                                                            \________/
%                                                                             1 2 3 4 5 6
%       T_Exp:                      n_t x n_p x n_{inv matrix containing the temperatures (^oC) used to constrain the inversion.
%                                         Each column of the matrix represents the n_t temperatures (in chronological order) corresponding to one of the n_p probes 
%                                         while each page of the matrix contains the temperatures used to do an inversion. 
%       Qdot:                         Total fluid heating power (W), a n_t x 1 x n_{inv variable used only to compute f for the FLSM. 
%                                         Each row corresponds to the mean heating power during an interval recorded at the end of each interval (in chronological order).
%       Vdot:                          Mean circulation flow rate (m^3/s) during the TRT (1 x 1 x n_{inv). 
%                                         The flow rate is assumed constant during the TRT and is used only to specify the fluid velocity in the vertical pipes of the TRCM.
%       delta_T:                     Temperature difference (^oC) between the inlet and outlet of the BHE (n_t x 1 x n_{inv). 
%                                         This variable is used only to compute f for the TRCM.  Each row contains the mean delta T during a time interval while each page contains a vector of data 
%
% Para: Structure array containing the following fields
% 
%      H:           Borehole Length (m)
%      r_b:        Borehole - Radius (m) 
%      r_i:         Pipe - Inner radius (m)
%      r_o:        Pipe - Outer radius (m)
%      k_Pipe:  Pipe - Thermal conductivity (W/m/K)
%      C_Pipe: Pipe - Volumetric heat capacity (J/mˆ3/K)
%      C_Fluid: Fluid - Volumetric capacity (J/mˆ3/K)
%      Rf:         Fluid - Convective resistance (mK/W)
%
% Inversion: Structure array containing the following fields
%
%      Xb: 6 × 2 matrix containing the boundary values for each parameter. The left column specifies the minimum value while the right column
%            indicates its maximum value. It is possible to fix the value of a parameter by putting the same value on all the elements of a given line. 
%            During an inversion, this value will be used by TRT-SInterp. Each line of Xb corresponds to a parameter. From line 1 to 6, these parameters are:
%                             1 - ks - Ground - Thermal conductivity (W/m/K)
%                             2 - Cs - Ground - Volumetric heat capacity (J/mˆ3/K)
%                             3 - kg - Grout - Thermal conductivity (W/m/K)
%                             4 - Cg - Grout - Volumetric heat capacity (J/mˆ3/K)
%                             5 - D - Half pipe spacing - Center to center (m)
%                             6 - T0 - Initial ground temperature (ˆoC)
%
%       FLSM and/or TRCM: By creating the fields FLSM and/or TRCM, the user can choose which interpretation model and which stopping criteria to use. 
%       These fields contain the following subfields:
%                                 SC_Inv : Specifies the stopping criteria of an inversion (^oC). 
%                                 SC Exp: Specifies the maximum relative change (%) allowed between two successive inversions of the Monte Carlo experiment. 
%
%% Output variables:
% X:    Contains the solution of each inverse problem and holds the subfields:
%  
%             Xb:       4 or 6  x 2 matrix containing the boundary on each parameter. For the subfield FLSM, lines 1 to 4 correspond to k_s,C_s,k_g and T_0 respectively 
%                          while for the field TRCM lines 1 to 6 correspond to k_s, C_s, k_g, C_g, D and T_0 respectively.
%             X0        4 or 6 x n_inv matrix containing the initial seed solution used by the inversion algorithm.  Each column corresponds to the random seed of an inversion.
%             Xs        4 or 6 x n_inv matrix containing the parameters found at the end of each inversion.  Notice that a given element of Xs corresponds to the same element of X0
%             Rb        2 x n_inv matrix of equivalent borehole resistances R_b (mK/W). The first row of the matrix is the borehole resistance obtained with the random thermal
%                         conductivity value of the grout (X0(3,:)).  The second row is the borehole resistance obtained with the optimized k_g value present in Xs(3,:).  
%                  
% T:     Holds the field Ts which represents the solution of temperature computed by the models using the values contained in Xs. 
%         For the FLSM, Ts is a n_t x 1 x n_inv matrix of the mean fluid temperatures while for TRCM Ts is a  n_t x n_p x n_inv matrix of temperatures at the 
%          location of the probes.
%
% Stats:  Contains the statistics of the inversions.
%
%                 mae             mean absolute error (^oC) between the measured and modeled temperatures.
%                 rmse            root-mean-square error (^oC), which corresponds to F.
%                 sse              sum of the squared error (^oC^2).
%                 CumMean   cumulative moving average of the parameters (same unit as Xb).
%                 CumStd       cumulative moving standard deviation of the parameters (same unit as Xb).
%                 CTime         computation time (s) for the inversion.
%                 nObjFun      number of direct problem solution (-).
% 
%
%% Sub-functions used
%
%  MAIN PROGRAM :
%       TRT_SInterp
%
% EXTERNAL FUNCTIONS
% 
%             BoreholeR
%             ConjGrad
%             FLS
%             Multipole_J
%             stop_TRT_ObjFunc
%             Tdot_TRCM_3D
%             TRCM_3D
%             TRCM_Spectral
%             TRT_Obj
%
%
%% KNOWN LIMITATIONS
%
%
%%%%%%%%%%%%%%%   MAIN PROGRAM  %%%%%%%%%%%%%%%%%
%
%% 1.0 - Initialization of variables

% 1.1 - Parameters of conjugate gradient/line search
Options.Tolerance_F=5e-3; % Tolerance of F
Options.TolRel_F=5e-3;       % Relative tolerance on F
Options.TolRel_t=5e-2;        % Relative tolerance of the step size alpha
Options.max_Iter=100;         % Maximum number of CG iterations

% 1.2 - Initialization of additional constants in Para
Para.L=0;
Para.r_Soil=10;

Para.rho_Soil=1000;
Para.rho_Grout=1000;
Para.Cp_Fluid=Para.C_Fluid/1000;
Para.rho_Fluid=1000;
Para.Cp_Pipe=Para.C_Pipe/1000;
Para.rho_Pipe=1000;

Para.Rp=log(Para.r_o/Para.r_i)/(2*pi*Para.k_Pipe)/2;
% Para.Vdot=Data.Vdot;

Para.Options.BC='TConstant';

% 1.3 - Preparation of the vector of times
Para.nt=size(Data.T_Exp,1);
Data.t=[Data.dt:Data.dt:Data.dt*Para.nt]';

% 1.4 - Preparation of the seed solution X0 for TRCM
n_inv=numel(Data.Vdot);
Xb=Inversion.Xb;

id=Xb(5,:)>Para.r_b-Para.r_o;   % Check if the pipes fit within the borehole
Xb(5,id)=Para.r_b-Para.r_o-1e-4;
id=Xb(5,:)<Para.r_o;                 % Check if there is an overlap with the pipes
Xb(5,id)=Para.r_o+1e-4;

TRCM.id=Xb(:,1)~=Xb(:,2);% Check if the lower and upper boundary are the same
TRCM.Xb=Xb;
TRCM.X0=repmat(mean(Xb,2),1,n_inv);

% Draw a random number from 0 to 1 and scale it to the permitted range of values.
TRCM.X0(TRCM.id,:)=repmat(Xb(TRCM.id,1),1,n_inv) + repmat(Xb(TRCM.id,2)-Xb(TRCM.id,1),1,n_inv).*rand(sum(TRCM.id),n_inv);

% 1.5 - Preparation of the seed solution X0 for FLSM - Neglect D and Cg.
FLSM.id=TRCM.id([1 2 3 6]);
FLSM.Xb=Xb([1 2 3 6],:);
FLSM.X0=TRCM.X0([1 2 3 6],:);

% 1.6 - Additional inputs
if isfield(Data,'id_F')
    Para.id_F=logical(Data.id_F);
else
    Para.id_F=true(numel(Data.t),1);
end

Data_0=Data;

Rb=@(Para,ks,kg,D)(1/(4*pi*kg)*(log(Para.r_b/Para.r_o)+log(Para.r_b/2./D)+(kg-ks)/(kg+ks)*log((Para.r_b./D).^4./((Para.r_b./D).^4-1)))+Para.Rp);

%% 2.0 - Stochastic Inversion with the finite line-source model

if isfield(Inversion,'FLSM')
    % 2.1 - Preparation of some variables
    Para.Options.Source='FLSM';
    Options.err_crit=Inversion.FLSM.SC_Inv;
    id=(Data.T_Probe_Location(2,:)==2 | Data.T_Probe_Location(2,:)==5);% Use only probes located in the fluid for FLSM
    
    % 2.2 - Merging new structs to Inversion
    Inversion.FLSM=cell2struct([struct2cell(Inversion.FLSM); struct2cell(FLSM)],[fieldnames(Inversion.FLSM); fieldnames(FLSM)],1);
    
    % 2.3 - Print header to screen
    fprintf('%9s  \n','--------------------------------------------FINITE LINE-SOURCE MODEL-------------------------------------------')
    fprintf('%9s  \n','                                       |            Realization            |    Cumulative moving average     |')
    fprintf('%9s','Inversion','NFunc','Time','MAE','ks','Cs','kg','Tg','ks','Cs','kg','Tg')
    fprintf('%9s \n','')
    
    i=1;
    
    while i<=n_inv
        tStart=tic; % Time counter for the current inversion
        Inversion.cpt=i;
               
        % 2.4 - Transform variable Data
        Data.t=Data_0.t;
        Data.T_Exp=Data_0.T_Exp(:,:,i);
        Data.Qdot=Data_0.Qdot(:,:,i);
        
        % 2.5 - Do inversion with FLSM
        [xs,Ts,G,nObjFun,nLSearch] =ConjGrad(@TRT_ObjFunc,Options,Para,Data,mean(Data.T_Exp(:,id),2),Inversion);
        
        % 2.6 - Save temperatures
        T.FLSM.Ts(:,:,i)=Ts;
        
        % 2.7 - Save Xs
        Xs=FLSM.X0(:,i);
        Xs(Inversion.FLSM.id)=xs;
        FLSM.Xs(:,i)=Xs;
        
        % 2.8 - Transform kg to Rb and save Xs 
        FLSM.Rb(:,i)=[Rb(Para,Inversion.FLSM.X0(1,i),Inversion.FLSM.X0(3,i),mean(Inversion.Xb(5,:),2)); Rb(Para,Xs(1),Xs(3),mean(Inversion.Xb(5,:),2))];
        
        % 2.9 - Compute stats
        Stats.FLSM.mae(i)=mean(abs(mean(Data.T_Exp(:,id),2)-T.FLSM.Ts(:,i)));
        Stats.FLSM.rmse(i)=sqrt(mean(mean((mean(Data.T_Exp(:,id),2)-T.FLSM.Ts(:,i)).^2)));
        Stats.FLSM.sse(i)=sum((mean(Data.T_Exp(:,id),2)-T.FLSM.Ts(:,i)).^2);
        Stats.FLSM.CumMean(:,i)=mean(FLSM.Xs,2);
        Stats.FLSM.CumStd(:,i)=std(FLSM.Xs,0,2);
        Stats.FLSM.CTime(i)=toc(tStart);
        Stats.FLSM.nObjFun(i)=nObjFun;
        Stats.FLSM.nLSearch(i)=nLSearch;
        Stats.FLSM.G(:,i)=G;
                
        % 2.10 -  Check for stopping criteria
        if Stats.FLSM.mae(i)<=Inversion.FLSM.SC_Inv || Inversion.FLSM.SC_Inv==0 % Go the next inversion
            S='%5d %10d %10.3g %8.2g %9.3g %10.3e %6.2f %8.2f %8.2f %10.3e %6.2f %8.2f  \n';
            fprintf(S,[i nObjFun Stats.FLSM.CTime(i) Stats.FLSM.mae(i) FLSM.Xs(:,i)' Stats.FLSM.CumMean(:,i)']);
            if i>1
                Stats.FLSM.RelChange(:,i)=(Stats.FLSM.CumMean(:,i)-Stats.FLSM.CumMean(:,i-1))./Stats.FLSM.CumMean(:,i)*100;
                % Check if parameters are sufficiently stable
                if max(abs(Stats.FLSM.RelChange(:,i)))<Inversion.FLSM.SC_Exp
                    break
                end
            end
            i=i+1;
        else % If CG converges in a local minimum or don't converge , draw a new seed and restart the inversion
            TRCM.X0(TRCM.id,i)=repmat(Xb(TRCM.id,1),1,1) + repmat(Xb(TRCM.id,2)-Xb(TRCM.id,1),1,1).*rand(sum(TRCM.id),1);
            FLSM.X0(:,i)=TRCM.X0([1 2 3 6],i);
            Inversion.FLSM.X0=FLSM.X0;
            disp('Convergence to a local minimum - Restarting the inversion')
        end
    end
    % 2.11 - Gather data to variable FLSM
    X.FLSM.Xb=Inversion.Xb;
    X.FLSM.X0=FLSM.X0;
    X.FLSM.X0(3,:)=TRCM.X0(3,:);
    X.FLSM.Xs=FLSM.Xs;
    X.FLSM.Rb=FLSM.Rb;
end

%% 3.0 - Stochastic Inversion with the thermal resistance and capacity model

if isfield(Inversion,'TRCM')
    % 3.1 - Preparation of some variables
    Para.Options.Source='TRCM';
    Options.err_crit=Inversion.TRCM.SC_Inv;
    
    % 3.2 - Merging new structs to Inversion
    Inversion.TRCM=cell2struct([struct2cell(Inversion.TRCM); struct2cell(TRCM)],[fieldnames(Inversion.TRCM); fieldnames(TRCM)],1);
    TRCM.X0=Inversion.TRCM.X0;
    
    % 3.3 - Print header to screen
    S='%5d %10d %9.3g %10.2g %8.3g %10.3e %6.3g %10.3e %6.2f %8.2f %8.2f %10.3e %6.2f %10.3e %6.2f %8.2f  \n';
    fprintf('%9s  \n','----------------------------------------------------------THERMAL RESISTANCE AND CAPACITY MODEL----------------------------------------------------')
    fprintf('%9s  \n','                                       |                      Realization                    |             Cumulative moving average              |')
    fprintf('%9s','Inversion','NFunc','Time','MAE','ks','Cs','kg','Cg','D','Tg','ks','Cs','kg','Cg','D','Tg')
    fprintf('%9s \n','')
    
    i=1;
    
    while i<=n_inv
        tStart=tic;% Time counter for the current inversion
        Inversion.cpt=i;
        
        % 3.4 - Transform variable Data
        Data.t=Data_0.t;
        Data.T_Exp=Data_0.T_Exp(:,:,i);
        Data.delta_T=Data_0.delta_T(:,:,i);
        Data.Vdot=Data_0.Vdot(:,:,i);       
        
        % 3.5 - Do inversion with TRCM
        [xs,Ts,G,nObjFun,nLSearch] =ConjGrad(@TRT_ObjFunc,Options,Para,Data,Data.T_Exp,Inversion);
        
        % 3.6 - Save temperatures
        T.TRCM.Ts(:,:,i)=Ts;
        
        % 3.7 - Save Xs and Rb
        Xs=TRCM.X0(:,i);
        Xs(Inversion.TRCM.id)=xs;
        TRCM.Xs(:,i)=Xs;
        TRCM.Rb(:,i)=[Rb(Para,TRCM.X0(1,i),TRCM.X0(3,i),TRCM.X0(5,i)); Rb(Para,Xs(1),Xs(3),Xs(5))];
        
        % 3.8 - Compute stats
        Stats.TRCM.mae(i)=mean(mean(abs(Data.T_Exp-T.TRCM.Ts(:,:,i))));
        Stats.TRCM.rmse(i)=sqrt(mean(mean((Data.T_Exp-T.TRCM.Ts(:,:,i)).^2)));
        Stats.TRCM.sse(i)=sum(sum(Data.T_Exp-T.TRCM.Ts(:,:,i)).^2);
        Stats.TRCM.CumMean(:,i)=mean(TRCM.Xs,2);
        Stats.TRCM.CumStd(:,i)=std(TRCM.Xs,0,2);
        Stats.TRCM.CTime(i)=toc(tStart);
        Stats.TRCM.nObjFun(i)=nObjFun;
        Stats.TRCM.nLSearch(i)=nLSearch;
        Stats.TRCM.G(:,i)=G;
        
        % 3.9 -  Check for stopping criteria
        if Stats.TRCM.mae(i)<=Inversion.TRCM.SC_Inv  || Inversion.TRCM.SC_Inv==0 % Go the next inversion
            fprintf(S,[i nObjFun Stats.TRCM.CTime(i) Stats.TRCM.mae(i) TRCM.Xs(:,i)' Stats.TRCM.CumMean(:,i)']);
            if i>1
                Stats.TRCM.RelChange(:,i)=(Stats.TRCM.CumMean(:,i)-Stats.TRCM.CumMean(:,i-1))./Stats.TRCM.CumMean(:,i)*100;
                % Check if parameters are sufficiently stable
                if max(abs(Stats.TRCM.RelChange(:,i)))<Inversion.TRCM.SC_Exp
                    break
                end
            end
            i=i+1;
        else % If CG converges in a local minimum, draw a new seed and restart the inversion
            TRCM.X0(TRCM.id,i)=repmat(Xb(TRCM.id,1),1,1) + repmat(Xb(TRCM.id,2)-Xb(TRCM.id,1),1,1).*rand(sum(TRCM.id),1);
            Inversion.TRCM.X0=TRCM.X0;
            disp('Convergence to a local minimum - Restarting the inversion')
        end
    end
    
    % 3.10 - Gather data to variable TRCM
    X.TRCM.Xb=Inversion.Xb;
    X.TRCM.X0=TRCM.X0;
    X.TRCM.Xs=TRCM.Xs;
    X.TRCM.Rb=TRCM.Rb;
end

save Inversion_TRT.mat X T Stats
end