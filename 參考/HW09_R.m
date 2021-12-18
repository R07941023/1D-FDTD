clear all;
clc;

%�Ŷ��ث׳]�w
xdim = 500; %X�bgrid���ث�
Ez = zeros(1,xdim); %�q���ث�
Hy = zeros(1,xdim); %�ϳ��ث�
x_source_position = 0.5*xdim;  %�i����m


%�ɶ��ث׳]�w
Steps = 500; 
c = 1; %���t
delta_x = 1;
delta_t = delta_x/c; %���F��K�p��A�Ndx,dt�P���t���O��1

%����P��L�ѼƳ]�w
epsilon_0 = 1; 
epsilon = epsilon_0 * ones(1,xdim); %���q�`��
mu_0 = 1; 
mu = mu_0 * ones(1,xdim); %�Ͼɲv


for n=1:1:Steps

	% Update Hy from Ez
    for i=1:1:xdim-1
        Hy(i)=(delta_t/(delta_x*mu(i)))*(Ez(i+1)-Ez(i)) + Hy(i) ;
    end
    
	% Hy Source
    Hy_g(n) =-exp(-((n-20)^2)/50);
    Hy(x_source_position) = Hy(x_source_position) + Hy_g(n);
        
    % Update Ez from Hy
    for i=2:1:xdim
        Ez(i)=(delta_t/(delta_x*epsilon(i)))*(Hy(i)-Hy(i-1)) + Ez(i);
    end
    
    % Ez Source
	Ez_g(n) = exp(-(n-20)^2/50);
	Ez(x_source_position) = Ez(x_source_position) +Ez_g(n);

    
    % plot
    subplot(2,1,1)
    plot(1:1:xdim,Ez,'color','r');
    titlestring=['Ez'];
    title(titlestring,'color','k');
    xlabel('x');
    ylabel('Ez');
    axis([0 xdim -2 2]);
 
    subplot(2,1,2)
	plot(1:1:xdim,Hy,'color','b');
    titlestring=['Hy'];
    title(titlestring,'color','k');
    xlabel('x');
    ylabel('Hy');
    axis([0 xdim -2 2]);

    getframe;
end