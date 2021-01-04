function [] = light_beam_deflection()
% light_beam_deflection : function to create an animation
% for the gavitation light beam deflection effect.
%
% Author & support nicolas.douillet (at) free.fr, 2007-2021.


%--- Computational parameters ---%
sz = 41; % size of the space grid, integer >= 3; default : 41
xc = 0;  % grid centre X coordinate; default : 0
yc = 0;  % grid centre Y coordinate; default : 0
space_shape_param = 0.04; % space shape parameter; default : 0.04

[nx,ny] = meshgrid(1:sz,1:sz);

u = 1:13; % animation time base vector; default : 1:13
u = [u fliplr(u(1,1:end-1))]; % loop back animation time vector
u = u(1,1:end-1);

%--- Display parameters ---%
star_sz_mult_coeff = 0.5; % star size multiplying coefficient; default : 0.5 
time_lapse = 0.25; % animation time lapse; default : 0.25
title_text = 'Light beams deflection as a function of gravitational field';
title_on = true; % enable / disable title
filename = 'light_beam_gravity_deflection.gif';
titleString = 'Gravity field intensity';
cmap = true; % enable / disable colormap
cbar = true; % enable / disable colorbar
az = -22; % azimutal view angle; default : -22
el = 80; % elevation view angle; default : 80

%--- Display settings ---%
h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);
axis tight manual;


for t = 1:length(u)
    
        
    sigma = space_shape_param*u(1,t);
    Phi = compute_space_shape(sz,xc,yc,-u(t),sigma); 
    
    mesh(Phi), hold on;         
    alpha 0;
    
    milieu = ceil(0.5*sz)*ones(2,1);
    
    W1 = constraint_parabole(-u(1,t),0,sqrt(0.25*((sz-1)^2 + (sz-1)^2)),0.25*pi);
    W2 = constraint_parabole(-u(1,t),0,sqrt(0.25*((sz-1)^2 + (sz-1)^2)),1.25*pi);
    
    %--- Compute interpolated altitude ---%
    W1z = griddata(ny-milieu(1,1),nx-milieu(2,1),Phi,W1(1,:),W1(2,:),'cubic');
    W2z = griddata(ny-milieu(1,1),nx-milieu(2,1),Phi,W2(1,:),W2(2,:),'cubic');        
    
    W1 = cat(1,W1,W1z);
    W2 = cat(1,W2,W2z);
    
    nb = 0.5*length(W1);    
    
    for s = 1:nb                
        
        if s > 1
            
            plot3([W1(1,s-1),W1(1,s)]+milieu(1,1),[W1(2,s-1),W1(2,s)]+milieu(2,1),[W1(3,s-1),W1(3,s)],'Color','y','Linewidth',2.5), hold on;
            plot3([W1(1,nb+s-1),W1(1,nb+s)]+milieu(1,1),[W1(2,nb+s-1),W1(2,nb+s)]+milieu(2,1),[W1(3,nb+s-1),W1(3,nb+s)],'Color','y','Linewidth',2.5), hold on;
            
            plot3([W2(1,s-1),W2(1,s)]+milieu(1,1),[W2(2,s-1),W2(2,s)]+milieu(2,1),[W2(3,s-1),W2(3,s)],'Color','y','Linewidth',2.5), hold on;
            plot3([W2(1,nb+s-1),W2(1,nb+s)]+milieu(1,1),[W2(2,nb+s-1),W2(2,nb+s)]+milieu(2,1),[W2(3,nb+s-1),W2(3,nb+s)],'Color','y','Linewidth',2.5), hold on;
            
        end
        
    end
    
    plot3(milieu(1,1),milieu(2,1),max(max(Phi)),'o','Color',[0 0 1],'MarkerSize',star_sz_mult_coeff*s,'Linewidth',star_sz_mult_coeff*s), hold on;
    
    ax = gca;
    ax.Clipping = 'off';
    set(ax,'Color',[0 0 0]);    
    
    axis([1 sz 1 sz -10 0]), axis off;
    
    if title_on
        title(title_text,'FontSize',20,'Color',[1 1 1],'Position',[4*sz/5,5*sz/4,0]);
    end
    
    if cmap
        colormap('hot');
        if cbar
            hcb = colorbar;
            set(hcb,'YDir','reverse');
            colorTitleHandle = get(hcb,'Title');
        end
        set(colorTitleHandle ,'String',titleString,'Color',[1 1 1],'Fontsize',16);
    else
        colormap([1 1 1]);
    end
    
    view(az,el);
    campan(0,2.3);
    
    drawnow;
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    %--- Write to the .gif file ---%
    if t == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',Inf,'DelayTime',time_lapse);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
    end
    
    clf;
    
end


end % light_beam_deflection


function [z] = compute_space_shape(sz, xc, yc, height, sigma)
%
% Author & support nicolas.douillet (at) free.fr, 2007-2021.


z = zeros(sz);

scale_min = -1;
scale_max = 1;
step = (scale_max-scale_min)/(sz-1);

xc = step*xc;
yc = step*yc;

sample_vect = scale_min:step:scale_max;
[i,j] = meshgrid(sample_vect);
idx_vect = round((0.5*(sz-1))*sample_vect+0.5*(sz-1)+1);

r = sqrt((i-xc).^2+(j-yc).^2);
z(idx_vect,idx_vect) = height*exp((-r.^2)/sigma^2);


end % compute_space_shape


function [W] = constraint_parabole(xmin, xmax, ymax, Alpha)
%
% Author & support nicolas.douillet (at) free.fr, 2007-2021.


step = 1e-1;
Xo = xmin:step:xmax;
r = ymax/sqrt(xmax-xmin);
c = -xmin;

Yo = r*sqrt(Xo+c);

W = [fliplr(Xo) Xo;fliplr(Yo) -Yo];

% Rotation matrix
Rm = @(Theta)[cos(Theta) sin(Theta);
              -sin(Theta) cos(Theta)];

W = Rm(Alpha)*W; % for x / y symetry


end % constraint_parabole