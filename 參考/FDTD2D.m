clear all
mu = pi*4e-7;
elp = 8.8e-12;
c = sqrt(1/mu/elp);
wavelength = 500e-9;
freq = c/wavelength;
dx = wavelength/20;
dy = dx;
dt = 4.5e-17;%1/c*((1/dx)^2+(1/dy)^2)^-0.5;

L = 300;
W = 300;

Ezx = zeros(W+1,L+1);
Ezy = zeros(W+1,L+1);
Bx = zeros(W+1,L+1);
By = zeros(W+1,L+1);

[x y] = meshgrid(1:L+1,1:W+1);
Ezx = 0.0005*exp(-((sqrt((x-L/2).^2+(y-W/2).^2)).^2)/(5^2));    %Guassian source
Ezy = 0.0005*exp(-((sqrt((x-L/2).^2+(y-W/2).^2)).^2)/(2^2));
% Bx = 0.5*exp(-((sqrt((x-L/2).^2)).^2)/(2^2))*elp/mu;
% By = -0.5*exp(-((sqrt((x-L/2).^2)).^2)/(2^2))*elp/mu;
Ez = Ezx + Ezy;

BBx = dt/mu/dx;
BBy = dt/mu/dy;
BEx = dt/elp/dx;
BEy = dt/elp/dy;

for timest = 1:1000
if rem(timest,10) == 1
    f = pcolor(Ez);
    colormap jet
    axis equal tight
    shading flat
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    xlabel('10 \mum');ylabel('10 \mum')
    caxis([-0.00005 0.00005])
%     rectangle('position',[(particle(1)-particle(3))/dx+pmln (particle(2)-particle(3))/dx 2*particle(3)/dx 2*particle(3)/dx],'curvature',1,'linewidth',2,'edgecolor',[1 0 0]);
    title(['t = ',num2str(timest*dt*10^15),' fs']);
    getframe;
end
%     hold on
%     rectangle('position',[W/2-50 L/2-50 101 101]);
%     hold off
%     axis([1 W+1 1 L+1]);
%     f = get(f,'cdata');
% %     colorbar
%     set(gca,'Units','normalized','position',[0.08 0.55 0.35 0.4])
%     title('with PML')
%     subplot(2,2,2)
%     imagesc(Evz(Wv/2-50:Wv/2+50,Lv/2-50:Lv/2+50));
%     axis([1 101 1 101]);
%     colorbar
%     caxis([min(min(f)) max(max(f))]);
%     caxis([-0.05 0.05]);
%     set(gca,'Units','normalized','position',[0.52 0.55 0.35 0.4])
%     title('without PML')
%     subplot(2,2,3)
%     plot(t)
%     title('RMSE')
%     ylabel('RMSE')
%     xlabel('time-step')
%     set(gca,'Units','normalized','position',[0.08 0.08 0.8 0.4])
%     set(gcf,'Units','normalized','position',[0.15 0.1 0.6 0.7])
%     getframe;
    
%     Ezx(:,1) = Ezx(:,1) + 0.5*sin(2*pi*freq*timest*dt);   %sin source
%     Ezy(:,1) = Ezy(:,1) + 0.5*sin(2*pi*freq*timest*dt);

    %update eq.
    Bx(1:W,1:L) = Bx(1:W,1:L) - BBy*(Ez(2:W+1,2:L+1) - Ez(2:W+1,1:L));
    By(1:W,1:L) = By(1:W,1:L) + BBx*(Ez(2:W+1,2:L+1) - Ez(1:W,2:L+1));
    Ezx(2:W+1,2:L+1) = Ezx(2:W+1,2:L+1) + BEx*(By(2:W+1,1:L) - By(1:W,1:L));
    Ezy(2:W+1,2:L+1) = Ezy(2:W+1,2:L+1) - BEy*(Bx(1:W,2:L+1) - Bx(1:W,1:L));
    Ez = Ezx + Ezy;
    
    %{
    if rem(timest,10) == 3
       h = gcf;
       fram = getframe(h);
       [im,map] = rgb2ind(fram.cdata,256,'nodither');
       if timest == 3
           imwrite(im,map,'PMLLECno.gif','WriteMode','overwrite','DelayTime',0.5,'LoopCount',inf);
       else
           imwrite(im,map,'PMLLECno.gif','WriteMode','append','DelayTime',0.5);
       end
    end
    if timest == 1000
        f = pcolor(Ez(pmln:W-pmln,1:L));
        colormap jet
        axis equal tight
        shading flat
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        xlabel('10 \mum');ylabel('10 \mum')
        rectangle('position',[(particle(1)-particle(3))/dx+pmln (particle(2)-particle(3))/dx 2*particle(3)/dx 2*particle(3)/dx],'curvature',1,'linewidth',2,'edgecolor',[1 0 0]);
        title(['t = ',num2str(timest*dt*10^15),' fs']);
        h = gcf;
        fram = getframe(h);
        [im,map] = rgb2ind(fram.cdata,256,'nodither');
        imwrite(im,map,'FDTD.gif','WriteMode','append','DelayTime',0.5);
    end
    %}
end