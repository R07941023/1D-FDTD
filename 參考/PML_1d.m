%Clearing variables in memory and Matlab command screen
clear all;
clc;
% Grid Dimension in x (xdim) direction
xdim=270;
%Total no of time steps
time_tot=350;
%Position of the source (center of the domain)
xsource=135;
%Courant stability factor
S=1;
% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;
% Spatial grid step length (spatial grid step = 1 micron)
delta=1e-6;
% Temporal gris step obtained using Courant condition
deltat=S*delta/c;
% Initialization of field vectors
Ez=zeros(1,xdim);
Hy=zeros(1,xdim);
% Initialization of permittivity and permeability vectors
epsilon=epsilon0*ones(1,xdim);
mu=mu0*ones(1,xdim);
%Perfectly matched layer boundary design:
% [Reference:-http://dougneubauer.com/wp-content/uploads/wdata/yee2dpml1/yee2d_c.txt
%  (An adaptation of 2-D FDTD TE code of Dr. Susan Hagness)]
%Value of boundary width in spacesteps
bound_width=35;
%Electrical conductivity vector initialization
sigma=zeros(1,xdim);
%Order of polynomial on which sigma is modeled
gradingorder=6;
%Required reflection co-efficient
refl_coeff=1e-6;
%Polynomial model for sigma
sigmamax=(-log10(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width*delta);
boundfact1=((epsilon(1,bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
boundfact2=((epsilon(1,xdim-bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
x=0:1:bound_width;
for i=1:1:xdim
    sigma(1,bound_width+1:-1:1)=boundfact1*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
    sigma(1,xdim-bound_width:1:xdim)=boundfact2*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
end
%Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
%for no reflection at every point
sigma_star=(sigma.*mu)./epsilon;
%Choice of natue of source
gaussian=1;
sine=0;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
impulse=0;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step
%Multiplication factor vectors for H vector update to avoid being calculated many times 
%in the time update loop so as to increase computation speed
A=((mu-0.5*deltat*sigma_star)./(mu+0.5*deltat*sigma_star)); 
B=(deltat/delta)./(mu+0.5*deltat*sigma_star);
                          
%Multiplication factor vectors for E vector update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                          
C=((epsilon-0.5*deltat*sigma)./(epsilon+0.5*deltat*sigma)); 
D=(deltat/delta)./(epsilon+0.5*deltat*sigma);                     
                   
% Update loop begins
for n=1:1:time_tot
    
    %if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ez(xsource)=1;
    end
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % vector where the wave has reached to avoid unnecessary updates.
    if n<xsource-2
        n1=xsource-n-1;
    else
        n1=1;
    end
    if n<xdim-1-xsource
        n2=xsource+n;
    else
        n2=xdim-1;
    end
    
    %Vector Update instead of for loop for Hy field incorporating magnetic
    %conductivity
    Hy(n1:n2)=A(n1:n2).*Hy(n1:n2)+B(n1:n2).*(Ez(n1+1:n2+1)-Ez(n1:n2));
    
    %Vector Update instead of for loop for Ez field incorporating magnetic
    %conductivity
    Ez(n1+1:n2+1)=C(n1+1:n2+1).*Ez(n1+1:n2+1)+D(n1+1:n2+1).*(Hy(n1+1:n2+1)-Hy(n1:n2));
    
    % Source conditions
    if impulse==0
        % If unit-time step
        if gaussian==0 && sine==0
            Ez(xsource)=1;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end
        %if gaussian
        if gaussian==1
            if n<=42
                Ez(xsource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
            end
        end
    else
        %if impulse
        Ez(xsource)=0;
    end
    %Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
    plot((1:1:xdim)*delta,Ez,'color','k');
    titlestring=['\fontsize{20}Plot of Ez vs x for 1D FDTD Perfectly Matched Layer boundary condition at time = ',num2str(round(n*deltat/10e-15)),' fs'];
    title(titlestring,'color','k');
    xlabel('x in m','FontSize',20);
    ylabel('Ez in V/m','FontSize',20);
    set(gca,'FontSize',20);
    axis([0 xdim*delta -3 3]);
    getframe;
end