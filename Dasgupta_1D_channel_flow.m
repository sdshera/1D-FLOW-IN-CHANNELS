%% Srijan Dasgupta - Assignment 3
clc
close all
clear all 
%% 1. Input Data 
% Physical Data 
Di=0.020;                                               %internal diameter (m)
Do=0.024;                                               %external diameter (m)
L=20;                                                   %length of tube (m)
x1=0;                                                   %start point of tube 
x2=L;                                                   %end point of tube
%Er(relative roughness)<0.0001 (smooth tube)
vin=1;                                                  %inlet velocity (m/s)
pin=2e5;                                                %inlet pressure (Pa) 
Tin=95+273.15;                                          %inlet temperature (K) 
roin=density(Tin,pin);                                  %inlet density 
Tw=95+273.15;                                           %wall temperature 
Text=20+273.15;                                         %ambient temperature 
pext=1e5;                                               %ambient pressure 
epsi=0.00001;                                           %relative roughness 
lambda_t=36;                                            %thermal conductivity (tube) carbon steel C 1.5%
% Numerical Data 
N=10;                                                   %number of control volumes
tolerance=1e-9;                                         %tolerance 
maxiter=1e6;                                            %maximum iteration 
%% 2. Previous Calculation
% Matrix definition 
x_cv=linspace(x1,x2,N+1);                               %positions of control volumes 
Xp_f=x_cv;                                              %positions of nodes (fluid) 
Xp_t=zeros(1,N+2);                                      %position matrix of nodes (tube)
Xp_t(1)=x1;                                             %1st node position       (m)
Xp_t(end)=x2;                                           %last node position      (m)
del_x=(x2-x1)/N;                                        %nodal distance (tube) 

for i=2:1:N+1
     Xp_t(i)=(x_cv(i)+x_cv(i-1))/2;                     %positions of nodes (tube) (m)
end 
% Geometry calculations 
S_f=(pi/4)*(Di^2);                                      %fluid flow cross-section 
S_t=(pi/4)*((Do^2)-(Di^2));                             %heat conduction cross-section 
Pi=pi*Di;                                               %inside perimeter 
Po=pi*Do;                                               %outside perimeter 
m_in=roin*vin*S_f;                                      %inlet mass flowrate (kg/s) (fluid) 

%% 3. Initial Temperature 
T_t=Text*ones(1,N+2);                                   %initial temperature of nodes (tube) 

%% 4a. Internal fluid Evaluation
v_f=vin*ones(1,N+1);                                    %inlet velocity (m/s)
P_f=pin*ones(1,N+1);                                    %inlet pressure (Pa) 
T_f=Tin*ones(1,N+1);                                    %inlet temperature (K) 
ro_f=roin*ones(1,N+1);                                  %inlet density 

v_star=v_f; 
P_star=P_f; 
T_star=T_f; 
ro_star=ro_f; 

error=ones(1,4);
T_prev=ones(1,N+2);                                     %Guess Temperature array
iter=1; 
residual=1;                                             %residual introduce 
while residual>tolerance && iter<maxiter 

%internal fluid loop
for i=1:1:N
    error=ones(1,4);
while max(error)>tolerance && iter<maxiter
    v_f(i+1)=v_f(i);                                    % estimated velocity at outlet 
    P_f(i+1)=P_f(i);                                    %estimated pressure at outlet 
    T_f(i+1)=T_f(i);                                    %estimated temperature at outlet 
    ro_f(i+1)=ro_f(i);                                  %estimated density at outlet 
    
    T_i(i)=(T_f(i+1)+T_f(i))/2;                         %control volume average temperature 
    P_i=(P_f(i+1)+P_f(i))/2;                            %control volume average pressure 
    v_i=(v_f(i+1)+v_f(i))/2;                            %control volume average velocity 
    ro_i=(ro_f(i+1)+ro_f(i))/2;                         %control volume average density 
    
    mu_i=viscosity(T_i(i),ro_i);                             %dynamic viscosity 
    mu_w_i=viscosity(T_t(i+2),ro_i); 
    cp_i=SP(T_i(i));                                    %specific heat 
    lambda_i=conductivity(T_i(i));                      %thermal conductivity 
    
    Re_i=(ro_i*v_i*Di)/mu_i;                            %reynolds number 
    Pr_i=(mu_i*cp_i)/lambda_i;                          %prandlt number 
    Gz_i=(Re_i*Pr_i*Di)/L;                              %graetz number
    alpha_i=alpha(Re_i,Pr_i,Gz_i,mu_i,mu_w_i,lambda_i,Di); %heat transfer coefficient 
    alpha_i=alpha_i*ones(1,N); 
    f_i=f(Re_i,epsi,Di);                                %friction factor 
    m_i=m_in;                                           %mass flow rate 
    
    
    % Solving N-S and state equations
    P_f(i+1)=(-m_i*(v_f(i+1)-v_f(i))+P_f(i)*S_f-f_i*Pi*del_x*((ro_i*v_i^2)/2))/S_f; %eq 2---->P[i+1]
    T_f(i+1)=(m_i*cp_i*T_f(i)-(m_in/2)*(v_f(i+1)^2-v_f(i)^2)+alpha_i(i)*Pi*del_x*...
        (T_t(i+1)-T_i(i)))/(m_i*cp_i);                  %eq 3--->T[i+1]
    ro_f(i+1)=density(T_f(i+1),P_f(i+1));               %eq 4---> ro[i+1]
    v_f(i+1)=m_i/(ro_f(i+1)*S_f);                       %eq 1--->v[i+1]
    
    
    error(1)=abs(T_star(i+1)-T_f(i+1));                 %temperature error 
    error(2)=abs(P_star(i+1)-P_f(i+1));                 %pressure error 
    error(3)=abs(ro_star(i+1)-ro_f(i+1));               %density error 
    error(4)=abs(v_star(i+1)-v_f(i+1));                 %velocity error
    
    v_star=v_f;
    P_star=P_f; 
    T_star=T_f; 
    ro_star=ro_f; 
    
    iter=iter+1;
end 

end 

%% 4b. External heat transfer coefficient 
g=9.81;
beta=1/Text;
ro_ext=pext/(287*Text); 
mu_ext=1.458*10^-6*Text^1.5/(Text+110.4); 
cp_ext=1031.5-0.210*Text+4.143e-4*Text^2;
lambda_ext=(2.728e-3+7.776e-5*Text);

% external heat transfer coefficient 
for i=2:N+1
    Gr_ext(i)=g*beta*ro_ext^2*(T_t(i)-Text)*Do^3/mu_ext^2;
    Pr_ext=mu_ext*cp_ext/lambda_ext;
    Ra_ext(i)=Gr_ext(i)*Pr_ext;
    if Ra_ext(i)<10^9
        C=0.47;
        n=1/4;
        k=1;
    else
        C=0.1;
        n=1/3;
        k=1;
    end
     Nu_ext(i)=C*Ra_ext(i)^n*k;
     alpha_ext(i)=(lambda_ext*Nu_ext(i))/Do;
end 
%% 5. Solve the tube
ap=ones(1,N+2);
bp=ones(1,N+2);
aw=ones(1,N+2);
ae=ones(1,N+2);
aw(1)=0; 
bp(1)=0;
ae(N+2)=0;
bp(N+2)=0; 
%Node i=2:N+1 (Central Nodes)
for i=2:1:N+1    
    aw(i)=(lambda_t*S_t)/(Xp_t(i)-Xp_t(i-1));                               %central nodes 
    ae(i)=(lambda_t*S_t)/(Xp_t(i+1)-Xp_t(i));                               %central nodes 
    ap(i)=aw(i)+ae(i)+alpha_i(i-1)*Pi*del_x+alpha_ext(i)*Po*del_x;          %central nodes
    bp(i)=alpha_i(i-1)*T_i(i-1)*Pi*del_x;                                   %central nodes 
end 

T_prev=ones(1,N+2);                                                         %Guess Temperature array


    for i=2:1:N+1
            T_t(i)=(aw(i)*T_t(i-1)+ae(i)*T_t(i+1)+bp(i))/ap(i);
    end
    T_t(1)=T_t(2); 
    T_t(end)=T_t(end-1);
%% 6. Is max|Tt[i]-Tt*[i]|<tolerance?? if no---->go to 4a
    residual=max(abs(T_t-T_prev));
    T_prev=T_t;  
    iter=iter+1;
    
end


%% 7a. Final calculations  
Re_i 
Pr_i
f_i
alpha_i=alpha_i(end)
alpha_o=alpha_ext(end)
vout=v_f(end)                                                               %final outlet velocity
Pout=P_f(end)                                                               %final outlet pressure 
Tout=T_f(end)-273.15                                                        %final outlet temperature

%Qw=m_i*cp_final*(T_f(end)-Tin); 

%% 7b. Print Results 
figure(1)
plot(Xp_t,T_t); 
xlabel('Length (m)'); 
ylabel('Temperature (K)');
hold on
plot(Xp_f,T_f);
legend('Wall Temperatue','Fluid Temperature')

%% Functions 
%density
function sum=density(T,P) 
    sum=847.2+1.298*T-(2.657e-3)*T^2; %water
%    sum=1164.45-0.4389*T-(3.21e-3)*T^4; %therminol 66
%    sum=P/(287*T); % air density 
end 

% viscosity water
function sum=viscosity(T,ro) 
    if T<353
       sum=0.9149-(1.2563e-2)*T+(6.9182e-5)*T^2-(1.9067e-7)*T^3+(2.6275e-10)*T^4-(1.4474e-13)*T^5; %water 
%       sum=ro*exp(-(-16.096+(586.38/(T-210.65)))); %therminol 66 
%       sum=((2.5393e-5)*sqrt(T/273.15))/(1+(122/T)); % air density
    else 
       sum=(3.7471e-2)-(3.5636e-4)*T+(1.3725e-6)*T^2-(2.6566e-9)*T^3+(2.5766e-12)*T^4-(1e-15)*T^5; %water
%       sum=ro*exp(-(-16.096+(586.38/(T-210.65)))); %therminol 66 
%       sum=((2.5393e-5)*sqrt(T/273.15))/(1+(122/T)); % air density
    end
end 

%thermal conductivity (page 26)
function sum=conductivity(T) 
    sum=-0.722+7.168*10^-3*T-9.137*10^-6*T^2; %simplified water
%    sum=0.116+(4.9e-5)*T-(1.5e-7)*T^2; %therminol 66
%    sum=2.728e-3+7.776e-5*T; % air density 
end 

%specific heat  
function sum=SP(T) 
    sum=5648.8-9.140*T+14.21e-3*T^2; %water
%    sum=658+2.82*T+(8.97e-4)*T^2; %therminol 66
%    sum=1031.5-0.21*T+4.143e-4*T^2; % air density 
end 

% K_i value 
function sum=alpha(Re,Pr,Gz,mui,muwi,lambda,D) 
    if Re<2000 && Gz>10 %short tube 
        K=((Di/del_x)^(1/3))*((mui/muwi)^(0.14));
        C=1.86;m=1/3;n=1/3; 
        Nu=C*(Re^m)*(Pr^n)*K; 
        sum=(Nu*lambda)/D;
    elseif Re<2000 && Gz<10 %long tube 
        K=1;
        C=3.66;m=0;n=0; 
        Nu=C*(Re^m)*(Pr^n)*K; 
        sum=(Nu*lambda)/D;
    elseif Re>2000 && (0.6<Pr<100) %highly viscous liquids 
        K=(mui/muwi)^(0.14); 
        C=0.027;m=0.8;n=0.33; 
        Nu=C*(Re^m)*(Pr^n)*K; 
        sum=(Nu*lambda)/D;
    else %gases 
        K=1; 
        C=0.023;m=0.8;n=0.04; 
        Nu=C*(Re^m)*(Pr^n)*K; 
        sum=(Nu*lambda)/D;
    end 
end 

% friction factor (Table B7) 
function sum=f(Re,epsi,Di) 
        if Re<2000
            sum=16/Re;
        else
            Er=epsi/Di;
            A=(2.457*log(1/((7/Re)^0.9+0.27*Er)))^16;
            B=(37530/Re)^16;
            sum=2*((8/Re)^12+1/(A+B)^(3/2))^(1/12);
        end 
end 

% Nusselt number (Table B1) 
function sum=Nu(Ra) 
        if Ra<=1e3
            sum=0; 
        elseif Ra>1e9
            C=0.1;
            n=1/3;
            K=1; 
            sum=C*Ra^n*K; 
        else 
            C=0.47;
            n=1/4;
            K=1; 
            sum=C*Ra^n*K; 
        end 
end

