function [P_final,Corr_out]=NRtracking2(varargin)
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

	Fmean=mean(mean(F));
	Ftemp=F-ones([subsize, subsize])*Fmean;
	dF_temp=sum(sum(Ftemp.^2));
	dF=sqrt(dF_temp);








	while flag==0
		dx=X-x0;
		dy=Y-y0;
		Funcval=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y)
		difference=1e-4;
		% Ftemp=F;
		dd1=S(Ftemp,dF,coef,P(1)+2*difference,P(2),P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		dd2=S(Ftemp,dF,coef,P(1),P(2)+2*difference,P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		dd3=S(Ftemp,dF,coef,P(1),P(2),P(3)+2*difference,P(4),P(5),P(6),subsize,dy,dx,X,Y);
		dd4=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4)+2*difference,P(5),P(6),subsize,dy,dx,X,Y);
		dd5=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5)+2*difference,P(6),subsize,dy,dx,X,Y);
		dd6=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5),P(6)+2*difference,subsize,dy,dx,X,Y);
		d1=S(Ftemp,dF,coef,P(1)+difference,P(2),P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		d2=S(Ftemp,dF,coef,P(1),P(2)+difference,P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		d3=S(Ftemp,dF,coef,P(1),P(2),P(3)+difference,P(4),P(5),P(6),subsize,dy,dx,X,Y);
		d4=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4)+difference,P(5),P(6),subsize,dy,dx,X,Y);
		d5=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5)+difference,P(6),subsize,dy,dx,X,Y);
		d6=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5),P(6)+difference,subsize,dy,dx,X,Y);
		c1=S(Ftemp,dF,coef,P(1)-difference,P(2),P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		c2=S(Ftemp,dF,coef,P(1),P(2)-difference,P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		c3=S(Ftemp,dF,coef,P(1),P(2),P(3)-difference,P(4),P(5),P(6),subsize,dy,dx,X,Y);
		c4=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4)-difference,P(5),P(6),subsize,dy,dx,X,Y);
		c5=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5)-difference,P(6),subsize,dy,dx,X,Y);
		c6=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5),P(6)-difference,subsize,dy,dx,X,Y);
		cc1=S(Ftemp,dF,coef,P(1)-2*difference,P(2),P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		cc2=S(Ftemp,dF,coef,P(1),P(2)-2*difference,P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y);
		cc3=S(Ftemp,dF,coef,P(1),P(2),P(3)-2*difference,P(4),P(5),P(6),subsize,dy,dx,X,Y);
		cc4=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4)-2*difference,P(5),P(6),subsize,dy,dx,X,Y);
		cc5=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5)-2*difference,P(6),subsize,dy,dx,X,Y);
		cc6=S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5),P(6)-2*difference,subsize,dy,dx,X,Y);
		der1=(-dd1+8*d1-8*c1+cc1)/(12*difference);
		der2=(-dd2+8*d2-8*c2+cc2)/(12*difference);
		der3=(-dd3+8*d3-8*c3+cc3)/(12*difference);
		der4=(-dd4+8*d4-8*c4+cc4)/(12*difference);
		der5=(-dd5+8*d5-8*c5+cc5)/(12*difference);
		der6=(-dd6+8*d6-8*c6+cc6)/(12*difference);
		% der2=(d2-c2)/(2*difference);
		% der3=(d3-c3)/(2*difference);
		% der4=(d4-c4)/(2*difference);
		% der5=(d5-c5)/(2*difference);
		% der6=(d6-c6)/(2*difference);
		if der1==0
			der1=inf;
		elseif der2==0
			der2=inf;
		elseif der3==0
			der3=inf;
		elseif der4==0
			der4=inf;
		elseif der5==0
			der5=inf;
		elseif der6==0
			der6=inf;
		end
		% [der1,errest,finaldelta] = derivest(@(P1) S(Ftemp,dF,coef,P1,P(2),P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y),P(1));
		% [der2,errest,finaldelta] = derivest(@(P2) S(Ftemp,dF,coef,P(1),P2,P(3),P(4),P(5),P(6),subsize,dy,dx,X,Y),P(2));
		% [der3,errest,finaldelta] = derivest(@(P3) S(Ftemp,dF,coef,P(1),P(2),P3,P(4),P(5),P(6),subsize,dy,dx,X,Y),P(3));
		% [der4,errest,finaldelta] = derivest(@(P4) S(Ftemp,dF,coef,P(1),P(2),P(3),P4,P(5),P(6),subsize,dy,dx,X,Y),P(4));
		% [der5,errest,finaldelta] = derivest(@(P5) S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P5,P(6),subsize,dy,dx,X,Y),P(5));
		% [der6,errest,finaldelta] = derivest(@(P6) S(Ftemp,dF,coef,P(1),P(2),P(3),P(4),P(5),P6,subsize,dy,dx,X,Y),P(6));
		J=[der1,der2,der3,der4,der5,der6]

		% for i=1:6
		% 	for j=1:6


		% 	end
		% end
		
		
		% for i=1:subsize
		% 	for j=1:subsize
		% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
		% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
		% 	end
		% end



		% dx=reshape(dx,[subsize*subsize,1]);
		% dy=reshape(dy,[subsize*subsize,1]);

		% [temp]=WarpFunc(dx,dy,P);

		% xpp=x0+temp(:,1);
		% ypp=y0+temp(:,2);
		% xp=reshape(xpp,[subsize,subsize]);
		% yp=reshape(ypp,[subsize,subsize]);
		% for i=1:subsize
		% 	for j=1:subsize
		% 		coef_y=floor(yp(i,j));
		% 		coef_x=floor(xp(i,j));
		% 		coef_in(i,j,:)=coef(coef_y,coef_x,:);
		% 	end
		% end
		% H=J'*J;
		% dP=-H\J';
		dP=Funcval./J
		% for i=1:6
		% 	if abs(dP(i))>5
		% 		dP(i)=5*dP(i)/abs(dP(i));
		% 	end
		% end
		P=P+dP
		if abs(Funcval)<0.0004
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

			flag=1;
			G_deformed=InterpFunc(xp,yp);
			meshcompare(F,G_deformed)
			P_final=P;
			Corr_out=sum(sum((F-G_deformed).^2));
		end
	end

end


function out=S(Ftemp,dF,coef,P1,P2,P3,P4,P5,P6,subsize,dy,dx,X,Y)
	xp=zeros([subsize, subsize]);
	yp=zeros([subsize, subsize]);
	G=zeros([subsize, subsize]);
	for i=1:subsize
		for j=1:subsize
			xp(i,j)=P1+P3.*dy(i,j)+dx(i,j).*(P2+1.0)+X(i,j);
			yp(i,j)=P4+P5.*dx(i,j)+dy(i,j).*(P6+1.0)+Y(i,j);
		end
	end

	for i=1:subsize
		% fprintf('interp %d \n', i);
		for j=1:subsize
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			% G(i,j)=[1, y_dec, y_dec^2, y_dec^3]*a*[1; x_dec; x_dec^2; x_dec^3];
		end
	end

	Gmean=mean(mean(G));
	Gtemp=G-ones([subsize, subsize])*Gmean;
	dG_temp=sum(sum(Gtemp.^2));
	dG=sqrt(dG_temp);
	out1=sum(sum((Ftemp./dF-Gtemp./dG).^2));
	out=1-out1/2;
	% out=sum(sum((Ftemp-G).^2));
end