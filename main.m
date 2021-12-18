clear all
clc

%% dimenstion parameter
xdim=150;
dx=1e-9;  % [m]
Steps=1000;
PML_w=20;
PML_n=6;
PML_R=1e-6;  % reflection coefficient

%% Souce
source=30;
intensity = 3;
wide = 1;
const = 50; 

%%
epsilon0=8.85e-12;
u0=1.2566e-6;
c=3e8;
dt=dx/c;
Ez=zeros(1,xdim);
Hy=zeros(1,xdim);
Ez_g=zeros(1,xdim);
Hy_g=zeros(1,xdim);
globel_Hy = zeros(Steps,xdim); 
globel_Ez = zeros(Steps,xdim);

%% three materials
epsilon=epsilon0*ones(1,xdim);
u=u0*ones(1,xdim);
a_meterial = [1,50]; 
b_meterial = [51,100]; 
c_meterial = [101,150];
epsilon(1,a_meterial(1):a_meterial(2))=1*epsilon0;
epsilon(1,b_meterial(1):b_meterial(2))=11.5*epsilon0;
epsilon(1,c_meterial(1):c_meterial(2))=10*epsilon0;
u(1,a_meterial(1):a_meterial(2))=1*u0;
u(1,b_meterial(1):b_meterial(2))=1*u0;
u(1,c_meterial(1):c_meterial(2))=1*u0;
record = [80,300]; %caculation reflection (simulaiton) [80,330]

%% PML
% Ez conductivity
PML_maxsigma=(-log10(PML_R)*(PML_n+1)*epsilon0*c)/(2*PML_w*dx);
PML_boundary_l=((epsilon(1,PML_w)/epsilon0)*PML_maxsigma)/((PML_w^PML_n)*(PML_n+1));
PML_boundary_r=((epsilon(1,xdim-PML_w)/epsilon0)*PML_maxsigma)/((PML_w^PML_n)*(PML_n+1));
sigma=zeros(1,xdim);  
x=0:PML_w;
for i=1:xdim
    sigma(1,PML_w+1:-1:1)=PML_boundary_l*((x+0.5*ones(1,PML_w+1)).^(PML_n+1)-(x-0.5*[0 ones(1,PML_w)]).^(PML_n+1));
    sigma(1,xdim-PML_w:xdim)=PML_boundary_r*((x+0.5*ones(1,PML_w+1)).^(PML_n+1)-(x-0.5*[0 ones(1,PML_w)]).^(PML_n+1));
end
% Hy conductivity
sigma_s=(sigma.*u)./epsilon;

% Hy coefficient 
A=((u-0.5*dt*sigma_s)./(u+0.5*dt*sigma_s)); 
B=(dt/dx)./(u+0.5*dt*sigma_s);
                          
% Ez coefficient                        
C=((epsilon-0.5*dt*sigma)./(epsilon+0.5*dt*sigma)); 
D=(dt/dx)./(epsilon+0.5*dt*sigma);      


%% 1D-FDTD

for time=1:1:Steps
    time;
    % time boundary
    if time < source-2  % left side
        xi=source-time-1;
    else
        xi=1;
    end
    if time < xdim-1-source % right side
        xf=source+time;
    else
        xf=xdim-1;
    end
    
    % Update Hy from Ez
    for i = xi:xf
        Hy(i)=A(i).*Hy(i)+B(i).*(Ez(i+1)-Ez(i));
    end
    
    % Hy source
    Hy_g(time) = intensity*exp(-((time-source)/wide)^2/const);
    
    % Update Ez from Hy
    for i = xi:xf
        Ez(i+1)=C(i+1).*Ez(i+1)+D(i+1).*(Hy(i+1)-Hy(i));
    end
    
    
    % Ez source
    Ez_g(time) = intensity*exp(-((time-source)/wide)^2/const);
    Ez(source) = Ez(source) + Ez_g(time);

    globel_Ez(time,:) = Ez(1,:);
    globel_Hy(time,:) = Hy(1,:);
    
    %Reflection
    totalE(1,time) = sum(abs(globel_Ez(time,:)));
    totalH(1,time) = sum(abs(globel_Hy(time,:)));
    if time == record(1)
        %caculate the simulation value
        rE1=sum(abs(globel_Ez(time,a_meterial(1):a_meterial(2))));
        RE1 = rE1/totalE(1,time);
        tE1=sum(abs(globel_Ez(time,b_meterial(1):b_meterial(2))));
        TE1 = tE1/totalE(1,time);
        rH1=sum(abs(globel_Hy(time,a_meterial(1):a_meterial(2))));
        RH1 = rH1/totalH(1,time);
        tH1=sum(abs(globel_Hy(time,b_meterial(1):b_meterial(2))));
        TH1 = tH1/totalH(1,time);
        %caculation the fresnel value
        n1 = sqrt(epsilon(1,a_meterial(1))*u(1,a_meterial(1)));
        n2 = sqrt(epsilon(1,b_meterial(1))*u(1,b_meterial(1)));
        f_R1 = ((n1-n2)/(n1+n2))^2;
        %error (fresnel-simulation)/fresnel
        error_RE1 = abs(f_R1-RE1)/f_R1*100;
        error_RH1 = abs(f_R1-RH1)/f_R1*100;
        string = [' '];
        disp(string)
        string = ['In first reflection'];
        disp(string)
        string = ['Ez Simulation value is ', num2str(RE1*100), '%, fresnel value is ', num2str(f_R1*100), '%, error is ', num2str(error_RE1), '%'];
        disp(string)
        % string = ['Hy Simulation value is ', num2str(RH1*100), '%, fresnel value is ', num2str(f_R1*100), '%, error is ', num2str(error_RH1), '%'];
        % disp(string)
    end
    if time == record(2)
        %caculate the simulation value
        rE2=sum(abs(globel_Ez(time,b_meterial(1):b_meterial(2))));
        RE2 = rE2/tE1;
        tE2=sum(abs(globel_Ez(time,c_meterial(1):c_meterial(2))));
        TE2 = tE2/tE1;
        rH2=sum(abs(globel_Hy(time,b_meterial(1):b_meterial(2))));
        RH2 = rH2/tH1;
        tH2=sum(abs(globel_Hy(time,c_meterial(1):c_meterial(2))));
        TH2 = tH2/tH1;
        %caculation the fresnel value
        n2 = sqrt(epsilon(1,b_meterial(1))*u(1,b_meterial(1)));
        n3 = sqrt(epsilon(1,c_meterial(1))*u(1,c_meterial(1)));
        f_R2 = ((n2-n3)/(n2+n3))^2;
        %error (fresnel-simulation)/fresnel
        error_RE2 = abs(f_R2-RE2)/f_R2*100;
        error_RH2 = abs(f_R2-RH2)/f_R2*100;
        string = [' '];
        disp(string)
        string = ['In second reflection'];
        disp(string)
        string = ['Ez Simulation value is ', num2str(RE2*100), '%, fresnel value is ', num2str(f_R2*100), '%, error is ', num2str(error_RE2), '%'];
        disp(string)
        % string = ['Hy Simulation value is ', num2str(RH2*100), '%, fresnel value is ', num2str(f_R2*100), '%, error is ', num2str(error_RH2), '%'];
        % disp(string)
    end
    
    % plot
    subplot(3,1,1)
    plot((1:xdim)*dx,Ez,'color','r');
    titlestring=['1D FDTD PML [Vacuum(0-0.5)/GaN(0-1)/Si(1-1.5)]'];
    title(titlestring,'color','k');
    xlabel('x [m]');
    ylabel('Ez [V/m]');
    axis([0 xdim*dx -2 2]);
    grid on

    subplot(3,1,2)
    plot((1:xdim)*dx,Hy,'color','b');
    xlabel('x [m]');
    ylabel('Hy [V/m]');
    axis([0 xdim*dx -3e-3 3e-3]);
    grid on

    subplot(3,1,3)
    plot((1:xdim)*dx,sigma)
    xlabel('x [m]');
    ylabel('sigma');
    axis([0 xdim*dx 0 900000]);
    getframe;
end