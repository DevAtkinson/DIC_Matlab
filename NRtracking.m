function [P_final,Corr_out]=NRtracking(varargin)
	for i=1:nargin/2
		switch varargin{i*2-1}
		case 'subset size'
			subsize=varargin{i*2};
		case 'subset position'
			subpos=varargin{i*2};
		case 'undeformed image'
			F_in=varargin{i*2};
		case 'deformed image'
			G_in=varargin{i*2};
		case 'guess'
			Pinitial=varargin{i*2};
		case 'coef'
			coef=varargin{i*2};%warp function
		case 'warp2'
			Warp2_func=varargin{i*2};% Jacobian
		case 'warp3'
			Warp3_func=varargin{i*2};%warp matrix
		end
	end
	G=G_in;
	F=F_in(subpos.coords(1):subpos.coords(3),subpos.coords(2):subpos.coords(4));
	[r,c]=size(F);
	XX=subpos.coords(2):subpos.coords(2)+subsize-1;										% x values over subset
	YY=subpos.coords(1):subpos.coords(1)+subsize-1;										% y values over subset
	X=repmat(XX,r,1);
	Y=repmat(YY',1,c);
	% x0=subsize/2+subpos(1);													% x value for subset centre
	% y0=subsize/2+subpos(2);													% y value for subset centre
	x0=subsize/2+subpos.coords(2);													% x value for subset centre
	y0=subsize/2+subpos.coords(1);													% y value for subset centre
	P=Pinitial;
	[r_g,c_g]=size(G);
	[Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
	InterpFunc=griddedInterpolant(Xmesh',Ymesh',G','cubic');
	flag=0;
	while flag==0
		dx=X-x0;
		dy=Y-y0;
		% size(X)
		% size(dx)
		dx=reshape(dx,[subsize*subsize,1]);
		dy=reshape(dy,[subsize*subsize,1]);
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
		xp=reshape(xpp,[subsize,subsize]);
		yp=reshape(ypp,[subsize,subsize]);
		for i=1:subsize
			for j=1:subsize
				coef_y=floor(yp(i,j));
				coef_x=floor(xp(i,j));
				coef_in(i,j,:)=coef(coef_y,coef_x,:);
			end
		end
		dP=-(Hessian10(F,coef,P')')\(Jacobian10(F,coef,P')')
		P=P+dP'
		if norm(dP)<0.00001
			meshcompare(F,G_deformed)
			flag=1;
			G_deformed=InterpFunc(xp,yp);
			P_final=P;
			Corr_out=sum(sum((F-G_deformed).^2));
		end
	end

end