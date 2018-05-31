%
% LAMBPLOT
% Make some convergence plots
%
% lambplot(emats, ematw, plot);
%
% emats: errors from SW4
% ematw: errors from WPP
% plot: 1,2,...: type of plot
%
function lambplot(emats, ematw, plot)

% get some simulation info
[p dist scpu wcpu msol]=lambtest();

nres=7;
ndist=10;
% line width
lw = 2;
% strarting indices
bres=2;
bdist=1;

clf;
hold on;
if plot == 1
% rel error as func of cpu time
%  d = 10;
  for d=1:ndist
    h=loglog(scpu(bres:nres),emats(bres:nres,d)/msol(d),"b",wcpu(bres:nres),ematw(bres:nres,d)/msol(d),"r");
    set(h,"linewidth",lw);
  end
elseif plot == 2
% rel error as func of points per wave length
%  d = 10;
  for d=1:ndist
    h=loglog(p(bres:nres),emats(bres:nres,d)/msol(d),"b",p(bres:nres),ematw(bres:nres,d)/msol(d),"r");
    set(h,"linewidth",lw);
  end
elseif plot == 3
% rel error as func of distance
%  r = 7
   for r=bres:7
     h=semilogy(dist,emats(r,:)./msol(1,:),"b",dist,ematw(r,:)./msol(1,:),"r");
     set(h,"linewidth",lw);
   end
elseif plot == 4
% abs error as func of distance
%  r = 7
   for r=bres:7
     h=semilogy(dist,emats(r,:),"b",dist,ematw(r,:),"r");
     set(h,"linewidth",lw);
   end
end

hold off;