%% This file can be used to initialize the variables required by TRT_SInterp
% ( ) refer to units
% [ ] refer to matrix dimension

clear
close all

s=RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(s);

%% 1 - Creation of variable Para %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Para.H=150;             % Depth of borehole (m)
Para.r_b=0.076;         % Borehole radius (m)
Para.r_i=0.017;         % Pipe inner radius (m)
Para.r_o=0.022;         % Pipe outer radius (m)
Para.C_Fluid=4.2e6;     % Fluid volumetric capacity (J/m^3/K)
Para.k_Pipe=0.40;       % Pipe thermal conductivity (W/mK)
Para.C_Pipe=1.9e6;      % Pipe volumetric capacity (J/m^3/K)
Para.Rf=1e-6;           % Fluid thermal resistance (mK/W)

%% 2 - Creation of variable Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data,~,~]=xlsread(['c:\TRT_SInterp_v1p1\','TRT_Dataset.xlsx'],'Dataset','A4:E21604');

n=100;
Data.dt=30;                                                                % Time (s)                     [nt x  1 x  1]
Data.T_Probe_Location=[[0 150 0];[2 2 5]];                                 % Location of the probe (m/-)  [2  x np x  1]
Data.T_Exp=repmat(data(:,2:4),[1 1 n]);                                    % Measured temperatures (^oC)  [nt x np x ninv]
Data.Qdot=repmat(data(:,5),[1 1 n]);                                       % Measured heating power (W)   [nt x  1 x ninv]
Data.Vdot=repmat(25/1000/60,[1 1 n]);                                      % Circulation flowrate (m^3/s) [1  x  1 x ninv]
Data.delta_T=Data.Qdot./repmat(Para.C_Fluid*Data.Vdot,[size(Data.T_Exp,1) 1 1]);          % Temperature difference (^oC) [nt x  1 x ninv]

%% 3 - Creation of variable Inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Inversion.Xb=[...
    1.5 3.5;...             % ks - Ground - Thermal conductivity (W/m/K)
    2.1e6 2.6e6;...         % Cs - Ground - Volumetric heat capacity (J/m^3/K)
    1.0 2.0;...             % kg - Grout  - Thermal conductivity (W/m/K)
    2.25e6 2.25e6;...       % Cg - Grout  - Volumetric heat capacity (J/m^3/K)
    0.05 0.05;...           % D  - Half pipe spacing - Center to center (m)
    10 10];                 % T0 - Initial ground temperature (^oC)

Inversion.FLSM.SC_Inv=0.125; % Stopping criteria for the inversion (^oC)
Inversion.FLSM.SC_Exp=0.00; % Stopping criteria for the experiment (%)

% % Uncomment the next two lines to use also the TRCM
% Inversion.TRCM.SC_Inv=0.065; % Stopping criteria for the inversion (^oC)
% Inversion.TRCM.SC_Exp=0.00;  % Stopping criteria for the experiment (%)

%% 4 - Run TRT_SInterp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,T,Stats]=TRT_SInterp(Data,Para,Inversion);

%% 5 - Prepare some figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Color=[
    1.0000      0              0
    1.0000      0.4063         0
    1.0000      0.8125         0
    0.7813      1.0000         0
    0.3750      1.0000         0
    0           1.0000    0.0313
    0           1.0000    0.4375
    0           1.0000    0.8438
    0           0.7460    1.0000
    0           0.3333    1.0000
    0           0              1
    ];

ind=1;
t=[Data.dt:Data.dt:Data.dt*size(Data.T_Exp,1)]/3600;
id=(Data.T_Probe_Location(2,:)==2 | Data.T_Probe_Location(2,:)==5);% Use only probes located in the fluid for FLSM

if isfield(T,'FLSM')
    %% Figure 1
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    plot(t,mean(Data.T_Exp(:,id,ind),2),'linestyle',':','color','b','linewidth',2)
    plot(t,T.FLSM.Ts(:,:,ind),'linestyle','-','color','r','linewidth',2)
    xlabel('Time (h)')
    ylabel('Temperature (^oC)')
    box on
    
    H=legend('Experimental','FLSM');
    
    % Figure 2
    nc=25;
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1,'Parent',gcf,'FontSize',15);
    hold on
    hist(X.FLSM.Xs(1,:),nc)
    box on
    xlabel('k_s (W/mK)')
    ylabel('Frequency (-)')
    axis([2.6 3.5 0 20])
    axis square
    
    subplot(1,3,2,'Parent',gcf,'FontSize',20);
    hold on
    hist(X.FLSM.Xs(2,:)/1e6,nc)
    box on
    xlabel('C_s (MJ/Km^3)')
    ylabel('Frequency (-)')
    axis([1.5 3.0 0 20])
    axis square
    
    subplot(1,3,3,'Parent',gcf,'FontSize',20);
    hold on
    hist(X.FLSM.Xs(3,:),nc)
    box on
    xlabel('k_g (W/mK)')
    ylabel('Frequency (-)')
    axis([1.2 3.2 0 20])
    axis square
end

if isfield(T,'TRCM')
    %% Figure 3
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    for i=1:size(T.TRCM.Ts(:,id,ind),2);
        plot(t,Data.T_Exp(:,i,ind),'linestyle',':','color',Color(i,:),'linewidth',1)
        plot(t,T.TRCM.Ts(:,i,ind),'linestyle','-','color',Color(i,:),'linewidth',2)
    end
    xlabel('Time (h)')
    ylabel('Temperature (^oC)')
    box on
    
    H=legend('Experimental - Supply Temperature',...
        'Experimental - \downarrow - z=150 m',...
        'Experimental - Return Temperature',...
        'TRCM - Supply Temperature',...
        'TRCM - \downarrow - z=150 m',...
        'TRCM - Return Temperature');
    set(H,'Location','NorthWest','fontsize',12);
    
    
    % Figure 4
    nc=25;
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1,'Parent',gcf,'FontSize',15);
    hold on
    hist(X.TRCM.Xs(1,:),nc)
    box on
    xlabel('k_s (W/mK)')
    ylabel('Frequency (-)')
    axis([2.6 3.5 0 20])
    axis square
    
    subplot(1,3,2,'Parent',gcf,'FontSize',20);
    hold on
    hist(X.TRCM.Xs(2,:)/1e6,nc)
    box on
    xlabel('C_s (MJ/Km^3)')
    ylabel('Frequency (-)')
    axis([1.5 3.0 0 20])
    axis square
    
    subplot(1,3,3,'Parent',gcf,'FontSize',20);
    hold on
    hist(X.TRCM.Xs(3,:),nc)
    box on
    xlabel('k_g (W/mK)')
    ylabel('Frequency (-)')
    axis([1.2 3.2 0 20])
    axis square
end



