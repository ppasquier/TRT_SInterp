function [T,q]=Multipole_J(keq,kg,rb,rp,Tf,xym,xyg,Rpn,J)
%% Description :
%
% Function Multipole_J is a Matlab implementation of the Fortran program published by Bennet et al. (1987).  
%
% Author : Philippe Pasquier (philippe.pasquier@polymtl.ca)
% Multipole_J Version 1.1.0 (03-Mar-2011 10:00)
% Compatible with Matlab 8.2.0.701 (R2013b)
%.
%% References : 
% 1) Bennet, J., Claesson, J. & Hellström, G., 1987. Multipole Method to Compute the Conductive Heat Flows to and between Pipes in a Composite
%     Cylinder, Note on Heat Transfer 3, Lund, Sweden: University of Lund, Department of Building Technology and Mathematical Physics.
%
%% Variable  :
% keq    : Ground thermal conductivity (W/mK)
% kg      : Grout thermal conductivity (W/mK)
% rb       : Borehole radius (m)
% rp      :  Outer radius of the pipe [1 x npipe] (m)
% Tf      :  Fluid temperature of each pipe [1 x npipe] (degC)
% xym   : Coordinates (x, y) of the center  of each pipe  [npipe x 2] (m)
% xyg    : Coordinates (x, y) of the evaluation points [npoint x 2] (m)
% Rpn   : Thermal resistance of each pipe (Rp+Rf)  [1 x npipe] (mK/W)
% J        : Order of the multipoles [scalar] (-) 
%
% T        : Temperature at point xyg [degC]
% q        : Heat flux irradiating from each pipe [1 x npipe] (W/m)
%
%% Syntax:
% [T,q]=Multipole_J(keq,kg,rb,[0.1 0.3],[25 35],[1 1;2 5],[xg yg],[0.1 0.05],10);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Initialization of some variables

N=numel(rp);
Lambda_b=kg;
Lambda=keq;
rc=20;
Bc=0;
Beta=Rpn*2*pi*Lambda_b;
Tc=0;
xm=xym(:,1)';
ym=xym(:,2)';
xg=xyg(:,1);
yg=xyg(:,2);

z=xm+1i*ym;
zc=conj(z);
rm=abs(z);

sigma=(Lambda_b-Lambda)/(Lambda_b+Lambda);

%% 2 - What follows is a transcription of the original code of Bennet et al. 1987

rmn=nan(N,N);
zmn=nan(N,N);
PILAM=1/(2*pi*Lambda);
PILAMB=1/(2*pi*Lambda_b);
ALBETC=Bc+log(rc/rb);

for i=1:N
    for j=1:N
        rmn(i,j)=abs(z(i)-z(j));
        zmn(i,j)=z(j).*zc(i);
    end
end

Ro_diag=1/(2*pi*Lambda_b)*(Beta+log(rb./rp)+sigma*log(rb^2./(rb^2-rm.^2)))+1/(2*pi*Lambda)*(log(rc/rb)+Bc);
Ro=1/(2*pi*Lambda_b)*(log(rb./rmn)+sigma*log(rb^2./abs(rb^2-zmn)))+1/(2*pi*Lambda)*(log(rc/rb)+Bc);
Ro(1:(N+1):end)=Ro_diag;
Ko=inv(Ro);

ZRC=zeros(N,J);
RPZ=zeros(N,J);
RPZMN=nan(N,N,J);

for m=1:N
    RBM=rb^2/(rb^2-abs(z(m))^2);
    RKo(m,m)=PILAMB*(log(rb/rp(m))+Beta(m)+sigma*log(RBM))+PILAM*ALBETC;
    if J>=1
        for k=1:J
            RPMZN(m,m,k)=(rp(m)*zc(m)*RBM/rb^2).^k;
            if abs(z(m))~=0
                ZRC(m,k)=(z(m)/rc).^k;
                RPZ(m,k)=(rp(m)/z(m)).^k;
            end
        end
    end
    for n=1:N
        if m~=n
            PMK=z(n)-z(m);
            RMN=abs(PMK);
            RBM=rb.^2/abs(rb^2-z(n)*zc(m));
            PRB=rp(m)*zc(n)/(rb^2-zc(n)*z(m));
            RKo(m,n)=PILAMB*(log(rb/RMN)+sigma*log(RBM))+PILAM*ALBETC;
            if J>=1
                for k=1:J
                    RPZMN(m,n,k)=(rp(m)/PMK).^k;
                    RPMZN(m,n,k)=PRB.^k;
                end
            end
        end
    end
end

%% INITIAL VALUES OF THE ENERGY FLOWS AND MULTIPOLES

for m=1:N
    qbeg(m)=0;
    for n=1:N
        qbeg(m)=qbeg(m)+Ko(m,n)*(Tf(n)-Tc);
    end
    q(m)=qbeg(m);
end

for m=1:N
    for k=1:J
        P(m,k)=0;
        Pc(k)=0;
    end
end

if J>0
    for ite=1:100
        
        EPSMAX=0;
        for m=1:N
            PMMAX=0;
            for k=1:J
                PMK=0;
                for n=1:N
                    PRB=1/(rb^2-conj(z(n))*z(m));
                    KFAK=1;
                    for jj=1:J
                        if n~=m
                            PMK=PMK+P(n,jj)*RPZMN(n,m,jj)*RPZMN(m,n,k)*KFAK;
                        end
                        JPEND=min(jj,k);
                        KFAK1=KFAK;
                        KFAK2=1;
                        for jprim=0:JPEND
                            JJPRIM=jj-jprim;
                            KJPRIM=k-jprim;
                            TERM=1;
                            if JJPRIM>=1
                                TERM=conj(RPMZN(n,m,JJPRIM));
                            end
                            if KJPRIM>=1
                                TERM=TERM*RPMZN(m,n,KJPRIM);
                            end
                            PMK=PMK+conj(P(n,jj))*sigma*TERM*(rp(m)*rp(n)*PRB).^(jprim)*KFAK1*KFAK2;
                            if jprim~=JPEND
                                KFAK1=KFAK1*KJPRIM/(k+jj-1-jprim);
                                KFAK2=KFAK2*JJPRIM/(jprim+1);
                            end
                        end
                        KFAK=KFAK*(k+jj)/jj;
                    end
                    if n~=m
                        PMK=PMK+q(n)*PILAMB*RPZMN(m,n,k)/k;
                    end
                    PMK=PMK+q(n)*PILAMB*sigma*RPMZN(m,n,k)/k;
                end
                KFAK=1;
                for jj=k:J
                    PMK=PMK+Pc(jj)*(1-sigma)*RPZ(m,k)*ZRC(m,jj)*KFAK;
                    KFAK=KFAK*(jj+1)/(jj+1-k);
                end
                PMK=conj(PMK)*(Beta(m)*k-1)/(Beta(m)*k+1);
                PMMAX=max(abs(PMK),PMMAX);
                if abs(PMK)>1e-7
                    EPSMAX=max(EPSMAX,abs(PMK-P(m,k))/PMMAX);
                end
                P2(m,k)=PMK;
            end
        end
        
        %% NEW MULTIPOLES ARE ASSIGNED
        for m=1:N
            for k=1:J
                P(m,k)=P2(m,k);
            end
        end
        
        PMMAX=0;
        for k=1:J
            for m=1:N
                KFAK=1;
                for jj=1:k
                    PMK=PMK+P(m,jj)*ZRC(m,k)*RPZ(m,jj)*KFAK;
                    KFAK=KFAK*(k-jj)/jj;
                end
                PMK=PMK+ZRC(m,k)*q(m)*PILAMB/k;
            end
            XX=(1-Bc*k)/(1+Bc*k);
            PMK=conj(PMK)*XX*(sigma+1)/(sigma*XX*(rb/rc).^(2*k)-1);
            PMMAX=max(abs(PMK),PMMAX);
            if abs(PMK)>1e-7
                EPSMAX=max(EPSMAX,abs(PMK-Pc(k)))/PMMAX;
            end
            Pc(k)=PMK;
        end
        
        %% CALCULATION OF NEW ENERGY FLOWS
        
        for m=1:N
            QQQ=0;
            for n=1:N
                for jj=1:N
                    if m~=n
                        QQQ=QQQ+real(RPZMN(n,m,jj)*P(n,jj));
                    end
                    QQQ=QQQ+sigma*real(P(n,jj)*RPMZN(n,m,jj));
                end
            end
            for jj=1:J
                QQQ=QQQ+(1-sigma)*real(Pc(jj)*ZRC(m,jj));
            end
            QM(m)=QQQ;
        end
        
        for m=1:N
            QQQQ=0;
            for n=1:N
                QQQQ=QM(n)*Ko(m,n)+QQQQ;
            end
            q(m)=qbeg(m)-QQQQ;
        end
        
         
        q_ITE(ite,:)=q;
    end
end

% Evaluation of the temperatures
ZPR=complex(xg,yg);
for i=1:numel(xg)
    PMK=0;
    if abs(ZPR(i))<rb
        for n=1:N
            if J>0
                PRB=1/(rb^2-z(n)*conj(ZPR(i)));
                for k=1:J
                    PMK=PMK+P(n,k)*((rp(n)/(ZPR(i)-z(n))).^k+sigma*(rp(n)*conj(ZPR(i))*PRB).^k);
                end
            end
            RON=abs(ZPR(i)-z(n));
            PMK=PMK+q(n)*(PILAM*ALBETC+PILAMB*(log(rb/RON)+sigma*log(rb^2*abs(PRB))));
        end
        if J>0
            for k=1:J
                PMK=PMK+(1-sigma)*Pc(k)*(ZPR(i)/rc).^k;
            end
        end
    elseif  abs(ZPR(i))>rb && abs(ZPR(i))<rc
        for n=1:N
            if J>0
                for k=1:J
                    PMK=PMK+P(n,k)*(1+sigma)*(rp(n)/(ZPR(i)-z(n))).^k;
                end
            end
            RON=abs(ZPR(i)-z(n));
            PMK=PMK+q(n)*(PILAM*(ALBETC+sigma*log(rb/abs(ZPR(i))))+PILAMB*(1+sigma)*log(rb/RON));
        end
        if J>0
            for k=1:J
                PMK=PMK+Pc(k)*((ZPR(i)/rc).^k-sigma*(rb.^2/(rc*conj(ZPR(i)))).^k);
            end
        end
    end
    
    T(i)=real(PMK)+Tc;
end

T=T';
id=T>1e10;
T(id)=nan;
[xg yg T];
end