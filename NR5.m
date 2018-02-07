function NR5
	%function to create the functions that will give the symbolic jacobian and hessian of the G values

	numP=2;
	% F=F_in(1:subsize,1:subsize);
	% F = sym('F_%d_%d', [subsize subsize]);
	coef = sym('A_%d', [1 16])
	% G = sym('G_%d_%d', [subsize subsize]);
	P = sym('P_%d', [numP,1])
	% subcord=sym('subcord_%d', [2,1]);
	% syms d	W=A*X;
	% [r_g,c_g]=size(G);
	% [Xmesh,Ymesh]=meshgrid(1:1:subsize,1:1:subsize);
	% Int = (Xmesh',Ymesh',G','cubic')


									% y values over subset
	% XX=1:subsize;										% x values over subset
	% YY=1:subsize;										% y values over subset
	% X=repmat(XX,subsize,1);
	% Y=repmat(YY',1,subsize);
	
	% x0=subsize/2;													% x value for subset centre
	% y0=subsize/2;	

	% dx=X-x0;
	% dy=Y-y0;

	tic
	syms dx dy X Y;
	% xp=P(1)+P(3).*dy+dx.*(P(2)+1.0)+X;
	% yp=P(4)+P(5).*dx+dy.*(P(6)+1.0)+Y;
	xp=P(1)+dx.*(1.0)+X;
	yp=P(2)+dy.*(1.0)+Y;

	
	a=reshape(coef,[4,4]);
	% x_dec=mod(xp,1);
	% y_dec=mod(yp,1);
	x_dec=xp-floor(xp)
	y_dec=yp-floor(yp);
	G=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];

	J=jacobian(G,P);
	H=hessian(G,P);

	matlabFunction(J,'File','JacobianStandard_p2','Optimize',true,'Vars',{coef,P,dx,dy,X,Y});
	matlabFunction(H,'File','HessianStandard_p2','Optimize',true,'Vars',{coef,P,dx,dy,X,Y});


end