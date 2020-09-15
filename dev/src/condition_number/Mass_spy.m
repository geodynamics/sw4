a = load('Mass.txt');
mass = reshape(a,13*13*3,13*13*3);
mass_diag = zeros(13*13*3,13*13*3);
for i = 1:13*13
    mass_diag((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3) = ...
        mass((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3);
end
figure(1)
spy(mass)
hold on
spy(mass_diag,'r')
set(gca,'fontsize',24)
figure(2)
spy(mass((7-1)*13*3+1:(8-1)*13*3,(7-1)*13*3+1:(8-1)*13*3))
hold on
spy(mass_diag((7-1)*13*3+1:(8-1)*13*3,(7-1)*13*3+1:(8-1)*13*3),'r')
set(gca,'fontsize',24)


% a = load('Mass_pre.txt');
% mass_pre = reshape(a,13*13*3,13*13*3);
% figure(3)
% spy(mass_pre)
% set(gca,'fontsize',12)
% figure(4)
% spy(mass_pre((7-1)*13*3+1:(8-1)*13*3,(7-1)*13*3+1:(8-1)*13*3))
% set(gca,'fontsize',12)

