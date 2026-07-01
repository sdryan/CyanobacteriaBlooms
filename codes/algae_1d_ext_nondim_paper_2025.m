clear
clc
clf

Track=0;

L = 4e10;  %Domain Length
M = 1e5;  %Number of discretized points on the domain
count = 0;


dx = L/M;  %Spatial Step
dt = .001; %(1.0/24.0);   % Time Step
Time = 80.0;   % Total Time of Sim

%c1 = 4;
%c2 = 1;
%Original Parameters
%kmax = 1.0; %3.0*.67*(24/22); %4; %.028; %Numerator parameter in M-M knietics of nitrogen consumption by Mike (and Ana consumption of Phos., assumed the same).
kmax = .028; %.67/day converts to .028/hr, max absorbtion rate of chemical by Mike
qmax = .007; %max absoprtion rate of chemical by Ana
l0 = 10; % 10 microns avergae of lengths mike and ana
t0 = 1/kmax; %.2; % 10 um/(4.9 um/s wind speed)
rhoN0 = 1874.4; %char conc of Nitrogen
rhoP0 = 84.4;  %char conc of Phos
rhoC = 3.76; % char conc of cyano
alpha1 = 155; %diffusion rate of cyano
alpha2 = 333.33; %diffusion rate of chem
kn = 530; %mug/L half life of chemical
qn = 530;  %mug/L half life of chemical
v = 0; %Background fluid flow in domain
k = 0; %Effect of fluid flow on bacteria; or chem concentrations
delta = 0.0025; %Estimate of decay rate of bacteria
gamma = .0025; %Estimate of decay rate of chem
eta1 = .0091; % Consumption Rate of N by Mike
zeta1 = eta1; % Consumption Rate of N by Ana
eta2 = 0*.228; %Consumption Rate of P by Mike
zeta2 = .01*.228+0*eta2; %Consumption Rate of P by Ana
% Removing Competition Term for Now
beta1 = 0.0*.0005; %0.0001;%.0001; %competition coefficient Mike
beta2 = 0.0*.0005; %0.0001; %.0001; %competition coefficient Ana
CFL = alpha1*dt/(dx*dx)
f0 = 300;


% Non-dimensional Parameters (Note: kmax = 1, kn disasspears in nondim
% eqns)
alpha1_nd = alpha1/(l0*l0*kmax);
alpha2_nd = alpha2/(l0*l0*kmax);
delta_nd = delta/kmax;
k_nd = k/(l0*t0); 
v_nd = v*t0/l0;
qmax_nd = qmax/kmax;
qn_nd = qn/kn;
gamma_nd = gamma/kmax;
eta1_nd = eta1*rhoC/kmax;
zeta1_nd = zeta1*rhoC/kmax;
eta2_nd = eta2*rhoC/kmax;
zeta2_nd = zeta2*rhoC/kmax;
beta1_nd = beta1*rhoC/kmax;
beta2_nd = beta2*rhoC/kmax;
f0_nd = f0/(kmax*kn);

% Initialize key quantities
u = zeros(M+1,2);
P = zeros(M+1,2);
N = zeros(M+1,2);
r = zeros(M+1,1);
w = zeros(M+1,2);


 x1 = 0:dx:L;
% Initial Conditions

    % Set interval for bump function (a1,b1) using the affine
    % transformation to (-1,1)=(c1,d1) in the bump function
    % f(t) = c1 + (d1-c1)/(b1-a1)*(t-a1)
    a1 = L/4.0;
    b1 = 3*L/8.0;
    
for i = 1:M+1
    u(i,1) = .1;
   
    P(i,1) =  .01 + 0*x1(i)*(L-x1(i)); %.1+(i-1)*.1;
    w(i,1) = .1;
   %x1 = 0:dx:L;
    N(i,1) =  .01 + 0*x1(i)*(L-x1(i));
    
    if (i < M/4)  % f(x) as a step function (later make it a bump smooth function)
        f(i,1) = 0;
    elseif (i <= 3*M/8 && i > M/4)
       % f(i,1) = 1;  % Step function
        f(i,1) = f0_nd*exp(-1.0/(1.0-(-1.0+(2.0)/(b1-a1)*(x1(i)-a1))^2));  %Nondim function scaling is t0/kn
       %  f(i,1) = 10.0*exp(-1.0/(1.0-(x1(i)-5*L/16)^2)); % Bump Function (Smooth -- nice for analysis)
         %Scaled by t_0/k_N
    else
        f(i,1) = 0;
    end
end

% Time Loop
for k = 1:Time/dt
    % Define r
       for i = 2:M
            r1(i) = N(i,1)/(1+N(i,1)); %  Old dim wayc1*(N(i,1))/(1+c2*(N(i,1)));
            r2(i) = P(i,1)/(1+P(i,1))*(qmax/kmax); %Old dim c4*P(i,1))/(1+c5*(P(i,1)));
            f1(i,1) = f(i,1);
            f2(i,1) = f(i,1);
       end
       if (k > 20000)
           for i = 1:M
               f1(i,1) = 0.0;
           end
       end
      if (k > 20000)
           for i = 1:M
               f2(i,1) = 0.0;
           end
       end
  
    % Evolution
    for i = 2:M
        u(i,2) = u(i,1) + dt*(alpha1_nd*(u(i+1,1)-2.0*u(i,1)+u(i-1,1))/(dx*dx)+r1(i)*u(i,1)*(1.0-u(i,1))-beta1*u(i,1)*w(i,1)-delta_nd*u(i,1)+k_nd*v*(u(i+1)-u(i-1))/(2*dx));
        w(i,2) = w(i,1) + dt*(alpha1_nd*(w(i+1,1)-2.0*w(i,1)+w(i-1,1))/(dx*dx)+r2(i)*w(i,1)*(1.0-w(i,1))-beta2*u(i,1)*w(i,1)-delta_nd*w(i,1)+k_nd*v*(w(i+1)-w(i-1))/(2*dx));
        N(i,2) = N(i,1) + dt*(alpha2_nd*(N(i+1,1)-2.0*N(i,1)+N(i-1,1))/(dx*dx)-gamma_nd*N(i,1)-eta1_nd*N(i,1)*u(i,1)-zeta1_nd*N(i,1)*w(i,1)+f1(i,1)+k_nd*v*(N(i+1)-N(i-1))/(2*dx)); 
        P(i,2) = P(i,1) + dt*(alpha2_nd*(P(i+1,1)-2.0*P(i,1)+P(i-1,1))/(dx*dx)-gamma_nd*P(i,1)-eta2_nd*P(i,1)*u(i,1)-zeta2_nd*P(i,1)*w(i,1)+f2(i,1)+k_nd*v*(P(i+1)-P(i-1))/(2*dx));
        
    end
    
      % Enforce the Boundary Conditions %No Flux
    u(1,2) = u(2,2);
    u(M+1,2) = u(M,2);
    w(1,2) = w(2,2);
    w(M+1,2) = w(M,2);
    
    P(1,2) = P(2,2);
    P(M+1,2) = P(M,2);  
    N(1,2) = N(2,2);
    N(M+1,2) = N(M,2);
    
    % Reset Time
    for i = 1:M+1
        u(i,1) = u(i,2);
        w(i,1) = w(i,2);
        P(i,1) = P(i,2);
        N(i,1) = N(i,2);
    end
    
    
    x1 = 0:dx:L;
 
    count = count+1;
    if count == 1000
        count = 0;
   % %hold on
   subplot(2,3,1)
   plot(x1*l0/1e9,u(:,2)*rhoC,'ro-') %green=Mike
   axis([0 L*l0/1e9 0 2*rhoC])
   title('Population of Mike')
   xlabel('Location (km)')
   ylabel('Population  (\mug/L)')

   subplot(2,3,2)
   plot(x1*l0/1e9,w(:,2)*rhoC,'bo-') %black=Ana
   axis([0 L*l0/1e9 0 2*rhoC])
   title('Population of Ana')
   xlabel('Location (km)')
   ylabel('Population  (\mug/L)')

   subplot(2,3,3)
   plot(x1*l0/1e9,(w(:,2)+u(:,2))*rhoC,'mo-') %black=Ana
   axis([0 L*l0/1e9 0 2*rhoC])
   title('Population of Bacteria')
   xlabel('Location (km)')
   ylabel('Population (\mug/L)')

   Time1 = k*dt;
   xt = [.45*L*l0/1e9, .55*L*l0/1e9, .77*L*l0/1e9];
    yt = [7-.1, 7-.1,7-.1];
    str={'t = ',num2str(Time1*(t0)/24,'%.2f'),'days'};
    text(xt,yt,str,'FontSize',16,'Color','blue')

   subplot(2,3,4)
   plot(x1*l0/1e9,N(:,2)*rhoN0,'go-')%blue=Nitrogen
   axis([0 L*l0/1e9 0 20000])
   title('Concentration of Nitrogen')
   xlabel('Location (km)')
   ylabel('Concentration (\mug/L)')

   subplot(2,3,5)
   plot(x1*l0/1e9,P(:,2)*rhoP0,'co-')%red=Phosphorus
   axis([0 L*l0/1e9 0 5000])
   title('Concentration of Phosphorus')
   xlabel('Location (km)')
   ylabel('Concentration (\mug/L)')

     subplot(2,3,6)
   plot(x1*l0/1e9,(P(:,2)*rhoP0+N(:,2)*rhoN0),'ko-')%red=Phosphorus
   axis([0 L*l0/1e9 0 25000])
   title('Concentration of Chemical')
   xlabel('Location (km)')
   ylabel('Concentration (\mug/L)')
   %hold off

    drawnow
    
%       Track = Track + 1;
%out=sprintf('Algae_mov_%3.3d.jpg',Track);
%filename=sprintf('Algae_mov_%3.3d.jpg',Track);
%saveas(gcf,filename)
    end
end

% Output Still Images for Video File

