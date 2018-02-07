function [P_final,Corr_out,iterations]=NRtracking3(varargin)
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
	% size(G)
	% meshcompare(F_in,G_in)
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
	count=0;

	dx=X-x0;
	dy=Y-y0;

	while flag==0
		count=count+1;
		
		% [G,coef_out]=Gvalues(coef,P,dy,dx,X,Y);
		% [J,H]=jacHes(F,G,coef_out,P,dx,dy,X,Y);
		% [G,J,H,Funcval]=getJac(coef,P,dy,dx,X,Y,F);
		[G,J,H,Funcval]=getJac5(coef,P,dy,dx,x0,y0,F);
		% GGGGG=G_in(subpos.coords(1):subpos.coords(3),subpos.coords(2):subpos.coords(4));
		% meshcompare(G,GGGGG)
		% G_interp=InterpFunc(xp,yp);
		% meshcompare(G,G_interp)
		% Funcval
		% J
		% H
		dP=H\J;

		% P=P+dP';
		P=P+dP';
		% Funcval=S2(F,Ftemp,dF,coef,P,subsize,dy,dx,X,Y);
		[G,J,H,Funcval]=getJac5(coef,P,dy,dx,x0,y0,F);
		% Funcval=CorVal3(F,G);
		if ((norm(dP)<0.0001)&&(Funcval<0.005))||(count>1000)||(Funcval<0.005)%||(Funcval<0.1) %abs(Funcval)<0.0004
			% dx=reshape(dx,[subsize*subsize,1]);
			% dy=reshape(dy,[subsize*subsize,1]);


			% x00=reshape(x0,[r*c,1]);
			% y00=reshape(y0,[r*c,1]);
			% determine the warped pixel points #4.15
			% might change
			% xp=x0+dx.*(1+P(2,converge))+P(3,converge).*dy+P(1,converge);
			% yp=y0+dy.*(1+P(6,converge))+P(5,converge).*dx+P(4,converge);
			% xp=x0+P(1,converge)+dx;
			% yp=y0+P(2,converge)+dy;


			% [temp]=WarpFunc(dx,dy,P);
			% % size(temp)
			% xpp=x0+temp(:,1);
			% ypp=y0+temp(:,2);
			% xp=reshape(xpp,[subsize,subsize]);
			% yp=reshape(ypp,[subsize,subsize]);


			% Fmean=mean(mean(F));
			% Gmean=mean(mean(G));
			% Mean_Corr=sum(sum((F-Fmean - (G-Gmean)).^2))
			% F2=sqrt(sum(sum(F.^2)));
			% G2=sqrt(sum(sum(G.^2)));
			% Norm_Corr=sum(sum((F./F2-G./G2).^2))
			% Actual_corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2))

			flag=1;
			% G_deformed=InterpFunc(xp,yp);
			% meshcompare(F,G)
			P_final=P;
			% Corr_out=sum(sum((F-G_deformed).^2));
			% [G,J,H,Funcval]=getJac6(coef,P,dy,dx,x0,y0,F);
			% meshcompare(F,G)
			% Corr_out=CorVal(coef,P,dy,dx,x0,y0,F);
			Corr_out=Funcval;
			iterations=count;
			% fprintf('Iterations: %d\n',count);
		end
	end

end


function out=S2(F,Ftemp,dF,coef,P,subsize,dy,dx,X,Y)
	xp=zeros([subsize, subsize]);
	yp=zeros([subsize, subsize]);
	G=zeros([subsize, subsize]);
	numP=max(size(P));
	for i=1:subsize
		for j=1:subsize
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
		end
	end

	for i=1:subsize
		% fprintf('interp %d \n', i);
		for j=1:subsize
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
			% G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			%% G(i,j)=[1, y_dec, y_dec^2, y_dec^3]*a*[1; x_dec; x_dec^2; x_dec^3];
		end
	end
	out=sum(sum((F-G).^2))/(sum(sum(F)));

	% Gmean=mean(mean(G));
	% Gtemp=G-ones([subsize, subsize])*Gmean;
	% dG_temp=sum(sum(Gtemp.^2));
	% dG=sqrt(dG_temp);
	% out=sum(sum((Ftemp./dF-Gtemp./dG).^2));


	% out=1-out1/2;
	% out=sum(sum((Ftemp-G).^2));
end

function out=CorVal(coef,P,dy,dx,X,Y,F)
	[r,c]=size(F);
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));

	F2=sqrt(sum(sum(F.^2)));
	G2=sqrt(sum(sum(G.^2)));

	% Corr=sum(sum((F./F2-G./G2).^2));
	out=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
end

function out=CorVal2(F,G)

	Fmean=mean(mean(F));
	Gmean=mean(mean(G));

	% F2=sqrt(sum(sum(F.^2)));
	% G2=sqrt(sum(sum(G.^2)));

	% Corr=sum(sum((F./F2-G./G2).^2));
	out=sum(sum((F-Fmean - (G-Gmean)).^2));
	% out=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
end

function out=CorVal3(F,G)

	% Fmean=mean(mean(F));
	% Gmean=mean(mean(G));

	F2=sqrt(sum(sum(F.^2)));
	G2=sqrt(sum(sum(G.^2)));

	% Corr=sum(sum((F./F2-G./G2).^2));
	out=sum(sum((F./F2-G./G2).^2));
	% out=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
end


function [G,J,H,Corr]=getJac(coef,P,dy,dx,X,Y,F)
	% custom correlation coeficient that works for newton raphson
	[r,c]=size(F);
	Jcoef=-2/(sum(sum(F.^2)));
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end
	numP=max(size(P));
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c
			% coef_out(i,j,:)=coef(floor(yp(i,j)),floor(xp(i,j)),:);
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
					
			% G(i,j)=[1, y_dec, y_dec^2, y_dec^3]*a*[1; x_dec; x_dec^2; x_dec^3];
			% Jacky=(JacobianValues(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianValues(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			% size(Jacky)
			% size(J)
			J=J+(F(i,j)-G(i,j)).*Jacky';
			H=H+Jacky'*Jacky;
		end
	end
	J=J.*Jcoef;
	H=-Jcoef*H;
	Corr=sum(sum((F-G).^2))/(sum(sum(F)));
	% Fmean=mean(mean(F))
	% Gmean=mean(mean(G))
	% out1=sum(sum((F-Fmean - (G-Gmean)).^2))
	% F2=sqrt(sum(sum(F.^2)));
	% G2=sqrt(sum(sum(G.^2)));
	% out2=sum(sum((F./F2-G./G2).^2))
	% actual_corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2))


	% meshcompare(G,F)
	% meshcompare(F,G)
	% max(max((G./F)))

	% figure
	% surf(G)
	% figure
	% surf(F)
end

function [G,J,H,Corr]=getJac2(coef,P,dy,dx,X,Y,F)
	% zero mean sum of squared difference
	[r,c]=size(F);
	Jcoef=-2/(sum(sum(F.^2)));
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	% for i=1:r
	% 	for j=1:c
	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
	% 	end
	% end
	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));

	
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c
			% coef_out(i,j,:)=coef(floor(yp(i,j)),floor(xp(i,j)),:);

			% G(i,j)=[1, y_dec, y_dec^2, y_dec^3]*a*[1; x_dec; x_dec^2; x_dec^3];
			% Jacky=(JacobianValues(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% % size(Jacky)
			% % size(J)
			% J=J+(F(i,j)-G(i,j)).*Jacky';
			% H=H+Jacky'*Jacky;

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			H=H+(-Jacky)'*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess);
		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=2*J;
	H=2*H;
	Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
end

function [G,J,H,Corr]=getJac3(coef,P,dy,dx,X,Y,F)
	% normalised sum of squared difference
	[r,c]=size(F);
	Jcoef=-2/(sum(sum(F.^2)));
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	% for i=1:r
	% 	for j=1:c
	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
	% 	end
	% end
	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	% Fmean=mean(mean(F));
	% Gmean=mean(mean(G));

	F2=sqrt(sum(sum(F.^2)));
	G2=sqrt(sum(sum(G.^2)));

	
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			% J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');
			J=J+(F(i,j)/F2-G(i,j)/G2).*(-Jacky'/G2);
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			% H=H+(-Jacky)'*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess);
			H=H+(-Jacky')*(-Jacky) + (F(i,j)*G2/F2-G(i,j)).*(-Hess); %if transpose jacky outside of brackets it doesn't work

		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=2*J;
	H=2/(G2*G2)*H;
	Corr=sum(sum((F./F2-G./G2).^2));
end

function [G,J,H,Corr]=getJac4(coef,P,dy,dx,X,Y,F)
	% zero mean sum of squared difference
	[r,c]=size(F);
	Jcoef=-2/(sum(sum(F.^2)));
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	% for i=1:r
	% 	for j=1:c
	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
	% 	end
	% end
	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));

	JJ=zeros([numP,1]);
	HH=zeros([numP,numP]);

	for i=1:r
		for j=1:c


			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			JJ=JJ+Jacky';
			HH=HH+Hess;
		end
	end
	JJ=JJ./(r*c);
	HH=HH./(r*c);
	
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c
			% coef_out(i,j,:)=coef(floor(yp(i,j)),floor(xp(i,j)),:);

			% G(i,j)=[1, y_dec, y_dec^2, y_dec^3]*a*[1; x_dec; x_dec^2; x_dec^3];
			% Jacky=(JacobianValues(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% % size(Jacky)
			% % size(J)
			% J=J+(F(i,j)-G(i,j)).*Jacky';
			% H=H+Jacky'*Jacky;

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));

			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky'+JJ);
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			H=H+(-Jacky+JJ)'*(-Jacky+JJ)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess+HH);
			% H=H+(-Jacky+JJ)'*(-Jacky+JJ);
		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=2*J;
	H=2*H;
	Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
end

function [G,J,H,Corr]=getJac5(coef,P,dy,dx,X,Y,F)
	% zero-mean sum of squared difference
	[r,c]=size(F);
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');
			% J=J+(F(i,j)/F2-G(i,j)/G2).*(-Jacky'/G2);
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			% H=H+(-Jacky)'*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess);
			H=H+(-Jacky')*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess); %if transpose jacky outside of brackets it doesn't work
			% H=H+(-Jacky)'*(-Jacky) + (F(i,j)*G2/F2-G(i,j)).*(-Hess);

		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=2*J;
	H=2*H;
	% Corr=sum(sum((F./F2-G./G2).^2));
	Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
end

function [G,J,H,Corr]=getJac6(coef,P,dy,dx,X,Y,F)
	% normalised sum of squared difference
	[r,c]=size(F);
	Jcoef=-2/(sum(sum(F.^2)));
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	% for i=1:r
	% 	for j=1:c
	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
	% 	end
	% end
	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	% Fmean=mean(mean(F));
	% Gmean=mean(mean(G));

	F2=sqrt(sum(sum(F.^2)));
	G2=sqrt(sum(sum(G.^2)));

	
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			% J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');
			J=J+(F(i,j)/F2-G(i,j)/G2).*((-Jacky')/G2);
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			% H=H+(-Jacky)'*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess);
			H=H+(-Jacky'./G2)*(-Jacky./G2) + (F(i,j)/F2-G(i,j)/G2).*(-Hess)/G2; %if transpose jacky outside of brackets it doesn't work

		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=2*J;
	% H=2/(G2*G2)*H;
	H=2*H;
	Corr=sum(sum((F./F2-G./G2).^2));
end

function [G,J,H,Corr]=getJac7(coef,P,dy,dx,X,Y,F)
	% zero-mean sum of squared difference (full differentiation)
	[r,c]=size(F);
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));
	Jmean=0;
	Hmean=0;
	for i=1:r
		for j=1:c
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Jmean=Jmean+Jacky';
			Hmean=Hmean+Hess;
		end
	end
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky'+Jmean/(r*c));
			% J=J+(F(i,j)/F2-G(i,j)/G2).*(-Jacky'/G2);
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			% H=H+(-Jacky)'*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess);
			H=H+(-Jacky'+Jmean/(r*c))*(-Jacky+Jmean'/(r*c))+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess+Hmean/(r*c)); %if transpose jacky outside of brackets it doesn't work
			% H=H+(-Jacky)'*(-Jacky) + (F(i,j)*G2/F2-G(i,j)).*(-Hess);

		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=2*J;
	H=2*H;
	% Corr=sum(sum((F./F2-G./G2).^2));
	Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
end

function [G,J,H,Corr]=getJac8(coef,P,dy,dx,X,Y,F)
	% normalised cross correlation
	[r,c]=size(F);
	% Jcoef=-2/(sum(sum(F.^2)));
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	numP=max(size(P));

	for i=1:r
		for j=1:c
			% xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			% yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
		end
	end

	% for i=1:r
	% 	for j=1:c
	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
	% 	end
	% end
	for i=1:r
		for j=1:c
			a=reshape(coef(floor(yp(i,j)),floor(xp(i,j)),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			if numP==6
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
			elseif numP==7
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
			elseif numP==8
				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
			end
		end
	end
	% Fmean=mean(mean(F));
	% Gmean=mean(mean(G));

	% F2=sqrt(sum(sum(F.^2)));
	% G2=sqrt(sum(sum(G.^2)));
	bottom=sqrt(sum(sum(F.^2))*sum(sum(G.^2)));

	
	
	J=zeros([numP,1]);
	H=zeros([numP,numP]);
	for i=1:r
		% fprintf('interp %d \n', i);
		for j=1:c

			% Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			% Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X(i,j),Y(i,j)));
			Jacky=(JacobianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));
			Hess=(HessianStandard(coef(floor(yp(i,j)),floor(xp(i,j)),:),P',dx(i,j),dy(i,j),X,Y));

			J=J+F(i,j)*Jacky';
			H=H+F(i,j)*Hess;

			% J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');
			% J=J+(F(i,j)/F2-G(i,j)/G2).*((-Jacky')/G2);
			% check1=(-Jacky)'*(-Jacky)
			% check2=(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess)
			% H=H+(-Jacky)'*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess);
			% H=H+(-Jacky'./G2)*(-Jacky./G2) + (F(i,j)/F2-G(i,j)/G2).*(-Hess)/G2; %if transpose jacky outside of brackets it doesn't work

		end
	end
	% J=J.*Jcoef;
	% H=-Jcoef*H;
	J=J/bottom;
	% H=2/(G2*G2)*H;
	H=H/bottom;
	Corr=2*(1-sum(sum((F*G).^2))/bottom);
end