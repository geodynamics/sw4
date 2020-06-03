% Plot 

%close all

% N6 = load('N6.txt');
% n1_c = N6(1);
% n2_c = N6(2);
% n3_c = N6(3);
% n1_f = N6(4);
% n2_f = N6(5);
% n3_f = N6(6);
% 
% X1c = load('X1c.txt'); X1c = reshape(X1c,n1_c,n2_c,n3_c); 
% X2c = load('X2c.txt'); X2c = reshape(X2c,n1_c,n2_c,n3_c);
% X3c = load('X3c.txt'); X3c = reshape(X3c,n1_c,n2_c,n3_c);
% 
% X1f = load('X1f.txt'); X1f = reshape(X1f,n1_f,n2_f,n3_f);
% X2f = load('X2f.txt'); X2f = reshape(X2f,n1_f,n2_f,n3_f);
% X3f = load('X3f.txt'); X3f = reshape(X3f,n1_f,n2_f,n3_f);
% 
% u1_c = load('uc1.txt'); u1_c = reshape(u1_c,n1_c,n2_c,n3_c);
% u1_f = load('uf1.txt'); u1_f = reshape(u1_f,n1_f,n2_f,n3_f);
% u2_c = load('uc2.txt'); u2_c = reshape(u2_c,n1_c,n2_c,n3_c);
% u2_f = load('uf2.txt'); u2_f = reshape(u2_f,n1_f,n2_f,n3_f);
% u3_c = load('uc3.txt'); u3_c = reshape(u3_c,n1_c,n2_c,n3_c);
% u3_f = load('uf3.txt'); u3_f = reshape(u3_f,n1_f,n2_f,n3_f);
% 
% 
% % this is used to generate 2d mesh (fixed one direction and calculate another two)
% for jj = 13
%     Xci = zeros(n1_c,n3_c);
%     Xcj = zeros(n1_c,n3_c);
%     u1cj = zeros(n1_c,n3_c);
%     for j = 1:n3_c
%         for i = 1:n1_c
%             Xci(i,j) = X1c(i,jj,j);
%             Xcj(i,j) = X3c(i,jj,j);
%             u1cj(i,j) = u2_c(i,jj,j);
%         end 
%     end
%     Xfi = zeros(n1_f,n3_f);
%     Xfj = zeros(n1_f,n3_f);
%     u1fj = zeros(n1_f,n3_f);
%     for j = 1:n3_f
% 	    for i = 1:n1_f
%             Xfi(i,j) = X1f(i,2*jj-1,j);
%             Xfj(i,j) = X3f(i,2*jj-1,j);
%             u1fj(i,j) = u2_f(i,2*jj-1,j);
% 	    end 
%     end
% 	% find min and max of u31
%     clow1 = min(min(u1cj));
%     clow2 = min(min(u1fj));
%     clow = 0.95*min(clow1,clow2);
%     chigh1 = max(max(u1cj));
%     chigh2 = max(max(u1fj));
%     chigh = 0.95*max(chigh1, chigh2);
%     clev = linspace(clow, chigh, 10); % these are the contour levels. You can change the number of levels (11)
%     figure(1), h1c=contourf(Xci,Xcj,u1cj,clev,'EdgeColor','none');
%     hold on
%     h1c1=contour(Xci,Xcj,u1cj,clev,'LineColor','k');
%     axis([0 2*pi 0 2*pi])
%     h1f = contourf(Xfi,Xfj,u1fj,clev,'EdgeColor','none');
%     h1f1 = contour(Xfi,Xfj,u1fj,clev,'LineColor','k'); 
%     caxis([clow, chigh]); % fix the colors
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('z')
%     title('u2')
% end 


load('times.txt')
load('energy_interior.txt')
%load('energy_bdry.txt')
%load('traction_continuity.txt')
a = 1;
times = times(a:end);
energy_interior = energy_interior(a:end);
%energy_bdry = energy_bdry(a:end);
figure(3)
plot(times,(energy_interior-energy_interior(1))/energy_interior(1))
axis([0 120 -5e-14 2e-14])
set(gca,'fontsize',24)
xlabel('t')
ylabel('(E(t)-E(0))/E(0)')
%figure(4)
%plot(times,energy_bdry-energy_bdry(1),'-o')
