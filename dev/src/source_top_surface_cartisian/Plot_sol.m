% Plot 

clear

N6 = load('N6.txt');
n1_c = N6(1);
n2_c = N6(2);
n3_c = N6(3);

X1c = load('X1c.txt'); 
X2c = load('X2c.txt'); 
X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);

u31_c = load('u3_1.txt'); u31_c = reshape(u31_c,n1_c,n2_c,n3_c);
u32_c = load('u3_2.txt'); u32_c = reshape(u32_c,n1_c,n2_c,n3_c);
u33_c = load('u3_3.txt'); u33_c = reshape(u33_c,n1_c,n2_c,n3_c);
u34_c = load('u3_4.txt'); u34_c = reshape(u34_c,n1_c,n2_c,n3_c);
u35_c = load('u3_5.txt'); u35_c = reshape(u35_c,n1_c,n2_c,n3_c);

% this is used to generate 2d mesh (fixed one direction and calculate another two)
for jj = 101
    Xci = zeros(n1_c,n3_c);
    Xcj = zeros(n1_c,n3_c);
    u31cj = zeros(n1_c,n3_c);
    u32cj = zeros(n1_c,n3_c);
    u33cj = zeros(n1_c,n3_c);
    u34cj = zeros(n1_c,n3_c); 
    u35cj = zeros(n1_c,n3_c);  
    for j = 1:n3_c
        for i = 1:n1_c
            Xci(i,j) = X1c(i);
            Xcj(i,j) = X3c(i,jj,j);
            u31cj(i,j) = u31_c(i,jj,j);
            u32cj(i,j) = u32_c(i,jj,j);
            u33cj(i,j) = u33_c(i,jj,j);
            u34cj(i,j) = u34_c(i,jj,j);
            u35cj(i,j) = u35_c(i,jj,j);
        end 
    end
	% find min and max of u31
    clow1 = min(min(u31cj));
    clow = 0.95*clow1;
    chigh1 = max(max(u31cj));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(1), h1c=contourf(Xci,Xcj,u31cj,clev,'EdgeColor','none');
    hold on
    h1c1=contour(Xci,Xcj,u31cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.1')
    
    clow1 = min(min(u32cj));
    clow = 0.95*clow1;
    chigh1 = max(max(u32cj));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(2), h2c=contourf(Xci,Xcj,u32cj,clev,'EdgeColor','none');
    hold on
    h2c1=contour(Xci,Xcj,u32cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.2')
    
    clow1 = min(min(u33cj));
    clow = 0.95*clow1;
    chigh1 = max(max(u33cj));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(3), h3c=contourf(Xci,Xcj,u33cj,clev,'EdgeColor','none');
    hold on
    h3c1=contour(Xci,Xcj,u33cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.3')
    
    clow1 = min(min(u34cj));
    clow = 0.95*clow1;
    chigh1 = max(max(u34cj));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(4), h4c=contourf(Xci,Xcj,u34cj,clev,'EdgeColor','none');
    hold on
    h4c1=contour(Xci,Xcj,u34cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.4')
    
    clow1 = min(min(u35cj));
    clow = 0.95*clow1;
    chigh1 = max(max(u35cj));
    chigh = 0.95*chigh1;
    clev = linspace(clow, chigh, 11); % these are the contour levels. You can change the number of levels (11)
    figure(5), h5c=contourf(Xci,Xcj,u35cj,clev,'EdgeColor','none');
    hold on
    h5c1=contour(Xci,Xcj,u35cj,clev,'LineColor','k');
    axis([0 2000 0 1000])
    caxis([clow, chigh]); % fix the colors
    colorbar
    axis equal
    xlabel('x')
    ylabel('z')
    title('u3,t=0.5')
    
    figure(jj)
    h1 = mesh(Xci,Xcj,Xci*0);view(2)
    set(h1,'EdgeColor','b')
    axis equal
    xlabel('x')
    ylabel('z')
    title('y = 1000')
    %print mesh.png
end 


