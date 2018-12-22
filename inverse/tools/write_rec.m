% write_rec
% Usage: write_rec()
% Input:
function write_rec(x0, x1, dx, y0, y1, dy, fname)
  if nargin < 7
    fname = 'rec-commands.in';
  end
  fp=fopen(fname,"w");
  nx = ceil((x1-x0)/dx);
  ny = ceil((y1-y0)/dy);
  xc = linspace(x0,x1,nx+1);
  yc=linspace(y0,y1,ny+1);

  for qy=1:ny+1
      for qx=1:nx+1
	  fprintf(fp,"rec file=sta-%d-%d x=%e  y=%e depth=0 nsew=1 usgsformat=1 sacformat=0\n", ...
		  qx-1, qy-1, xc[qx], yc[qy]);
      end
  end
  fclose(fp);
end
