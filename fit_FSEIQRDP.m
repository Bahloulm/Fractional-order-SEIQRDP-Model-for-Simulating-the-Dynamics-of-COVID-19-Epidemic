function [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,q1] = fit_FSEIQRDP(Q,R,D,Npop,E0,I0,time,guess,varargin)
%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-4);  
p.addOptional('tolFun',1e-4); 
p.addOptional('Display','iter'); 
p.addOptional('dt',0.1);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display  = p.Results.Display ;
dt  = p.Results.dt ;
%% Options for lsqcurvfit
options=optimset('TolX',tolX,'TolFun',tolFun,'MaxFunEvals',800,'Display',Display);
%% Fitting the data
input = [Q;R;D];
time = time';
fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; 
t = tTarget(1):dt:tTarget(end);
modelFun1 = @FSEIQRDP_fitting; 

% call Lsqcurvefit
[Coeff] = lsqcurvefit(@(para,t) modelFun1(para,t),guess,tTarget(:)',input,zeros(1,numel(guess)),[1 3 1 1 2 3 2 2],options);
%% Write the fitted coeff in the outputs
alpha1 = abs(Coeff(1));
beta1 = abs(Coeff(2));
gamma1 = abs(Coeff(3));
delta1 = abs(Coeff(4));
Lambda1 = abs(Coeff(5:6));
Kappa1 = abs(Coeff(7:8));
q1 = abs(Coeff(9:15));

 function [output] = FSEIQRDP_fitting(para,t0)
        alpha = abs(para(1));
        beta = abs(para(2));
        gamma = abs(para(3));
        delta = abs(para(4));
        lambda0 = abs(para(5:6));
        lambda = lambda0(1)*(1-exp(-lambda0(2).*t));
        kappa0 = abs(para(7:8));
        kappa = kappa0(1)*exp(-kappa0(2).*t);
        q=abs(para(9:15));
%% Initial conditions
N = numel(t);
Ef(1)=E0;
If(1)=I0;
Qf(1)=Q(1);
Rf(1)=R(1);
Sf(1)= Npop-Q(1)-R(1)-D(1)-E0-I0;
Df(1)=D(1);
Pf(1)=0;
qs=q(1);cps=1; 
qe=q(2);cpe=1;
qi=q(3);cpi=1;
qq=q(4);cpq=1;
qr=q(5);cpr=1;
qd=q(6);cpd=1;
qp=1;cpp=1;
dt = median(diff(t));
        for j=1:N
            c1(j)=(1-(1+qs)/j)*cps;
            cps=c1(j);

            c2(j)=(1-(1+qe)/j)*cpe;
            cpe=c2(j);

            c3(j)=(1-(1+qi)/j)*cpi;
            cpi=c3(j);

            c4(j)=(1-(1+qq)/j)*cpq;
            cpq=c4(j);

            c5(j)=(1-(1+qr)/j)*cpr;
            cpr=c5(j);

            c6(j)=(1-(1+qd)/j)*cpd;
            cpd=c6(j);

            c7(j)=(1-(1+qp)/j)*cpp;
            cpp=c7(j);
        end

        for i=2:N
             Sf(i) = dt^qs*( -alpha^qs*Sf(i-1) - (beta^qs./Npop).*If(i-1)*Sf(i-1)) -memo(Sf,c1,i) ;
             Ef(i) = dt^qe*( -gamma^qe*Ef(i-1) + (beta^qe./Npop).*If(i-1)*Sf(i-1)) -memo(Ef,c2,i) ;
             If(i) = dt^qi*( gamma^qi*Ef(i-1) - delta^qi.*If(i-1)) -memo(If,c3,i) ;
             Qf(i) = dt^qq*(  delta^qq.*If(i-1) -lambda(i-1)^qq*Qf(i-1) - kappa(i-1)^qq.*Qf(i-1)) -memo(Qf,c4,i) ;
             Rf(i)= (dt^qr)*(lambda(i-1)^qr)*Qf(i-1)-memo(Rf,c5,i);
             Df(i)= (dt^qd)*(kappa(i-1)^qd)*Qf(i-1)-memo(Df,c6,i);
             Pf(i)= (dt^qp)*(alpha^qp)*Sf(i-1)-memo(Pf,c7,i);
        end
    Q1 = Qf;
    R1 = Rf;
    D1 = Df;
    Q1 = interp1(t,Q1,t0);
    R1 = interp1(t,R1,t0);
    D1 = interp1(t,D1,t0);
    output = [Q1;R1;D1];
 end
end
