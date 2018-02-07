function out=lsqDIC(F_in,G,subpos,subsize,P)
	F=F_in(subpos.coords(1):subpos.coords(3),subpos.coords(2):subpos.coords(4));
	[r,c]=size(F);
	XX=subpos.coords(2):subpos.coords(2)+subsize-1;										% x values over subset
	YY=subpos.coords(1):subpos.coords(1)+subsize-1;										% y values over subset
	X=repmat(XX,r,1);
	Y=repmat(YY',1,c);
	% x0=subsize/2+subpos(1);													% x value for subset centre
	% y0=subsize/2+subpos(2);													% y value for subset centre
	x0=subsize/2+subpos.coords(2);													% x value for subset centre
	y0=subsize/2+subpos.coords(1);	
	[r_g,c_g]=size(G);
	[Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
	InterpFunc=griddedInterpolant(Xmesh',Ymesh',G','cubic');
	dx=X-x0;
	dy=Y-y0;
	% size(X)
	% size(dx)
	dx=reshape(dx,[r*c,1]);
	dy=reshape(dy,[r*c,1]);
	% x00=reshape(x0,[r*c,1]);
	% y00=reshape(y0,[r*c,1]);
	% determine the warped pixel points #4.15
	% might change
	% xp=x0+dx.*(1+P(2,converge))+P(3,converge).*dy+P(1,converge);
	% yp=y0+dy.*(1+P(6,converge))+P(5,converge).*dx+P(4,converge);
	% xp=x0+P(1,converge)+dx;
	% yp=y0+P(2,converge)+dy;
	[temp]=WarpFunc(dx,dy,P);
	% size(temp)
	xpp=x0+temp(:,1);
	ypp=y0+temp(:,2);
	xp=reshape(xpp,[r,c]);
	yp=reshape(ypp,[r,c]);

	% interpolate the investigated subset to obtain the subset used for comparison purposes
	G_deformed=InterpFunc(xp,yp);
	out=F-G_deformed;
end