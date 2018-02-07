function [P_final,Corr_out]=DICtracking2(varargin)
	%handle the input variables
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
		case 'warp'
			Warp_func=varargin{i*2};%warp function
		case 'warp2'
			Warp2_func=varargin{i*2};% Jacobian
		case 'warp3'
			Warp3_func=varargin{i*2};%warp matrix
		case 'correlation'
			corrleation_type=varargin{i*2};%what type of correlation to do
		end
	end

	G=G_in;
	%define reference subset
	% F=F_in(subpos(2):subpos(2)+subsize-1,subpos(1):subpos(1)+subsize-1);
	F=F_in(subpos.coords(1):subpos.coords(3),subpos.coords(2):subpos.coords(4));
	[r,c]=size(F);
	%assign subset from reference image for interpolation
	F_interp=F_in(subpos.coords(1)-1:subpos.coords(1)+subsize,subpos.coords(2)-1:subpos.coords(2)+subsize);	

	% XX=subpos(1):subpos(1)+subsize-1;										% x values over subset
	% YY=subpos(2):subpos(2)+subsize-1;										% y values over subset
	XX=subpos.coords(2):subpos.coords(2)+subsize-1;										% x values over subset
	YY=subpos.coords(1):subpos.coords(1)+subsize-1;										% y values over subset
	X=repmat(XX,r,1);
	Y=repmat(YY',1,c);
	% x0=subsize/2+subpos(1);													% x value for subset centre
	% y0=subsize/2+subpos(2);													% y value for subset centre
	x0=subsize/2+subpos.coords(2);													% x value for subset centre
	y0=subsize/2+subpos.coords(1);													% y value for subset centre
	[r_int,c_int]=size(F_interp);
	
	sumF=0;
	sumG2=0;
	Fmean=mean(mean(F)); % #4.3
	H=0;
	dF_temp=0;
	% determine the gradients of the images
	[gradx,grady]=imgradientxy(F_interp,'prewitt');

	xpmax=0;
	ypmax=0;
	% determine the hessian matrix
	for i=1:r
		for j=1:c
			%xp and yp are the distance of the current pixel point from the subset centre
			xp=XX(j)-x0;
			yp=YY(i)-y0;
			if xp>xpmax
				xpmax=xp;
			end
			if yp>ypmax
				ypmax=yp;
			end
			
			dfdx=gradx(i+1,j+1);
			dfdy=grady(i+1,j+1);
			%change
			% dWdp=[1 xp yp 0 0 0; 0 0 0 1 xp yp]; % Jacobian
			% dWdp=[1 0; 0 1];
			dWdp=dPFunc(xp,yp);

			dfdw{i,j}=[dfdx, dfdy]*dWdp;
			% dfdw_vec()=dfdw{i,j}
			H=H+dfdw{i,j}'*dfdw{i,j}; % #4.21
			dF_temp=dF_temp+(F(i,j)-Fmean)^2;
		end
	end
	dF=sqrt(dF_temp);
	P=Pinitial;

	[r_g,c_g]=size(G);
	[Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
	[r,c]=size(F);
	% Xmesh2=permute(Xmesh,[2,1]);
	% Ymesh2=permute(Ymesh,[2,1]);
	% G2=permute(G,[2,1]);
	% InterpFunc=griddedInterpolant(Xmesh2,Ymesh2,G2,'cubic');
	InterpFunc=griddedInterpolant(Xmesh',Ymesh',G','cubic');

	converge=1;
	flag=1;
	% perform the subset matching iteratively
	while (flag==1)
		% tic
		% distance of the current pixel point from the subset centre
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
		G_defromed=InterpFunc(xp,yp);
		% G_defromed=interp2(Xmesh,Ymesh,G,xp,yp,'cubic');

		% G_store{converge}=G_defromed;
		% meshcompare(F,G_defromed)
		G_def_mean=mean(mean(G_defromed)); % #4.4

		

		% determine the summed part of equation #4.20
		q=0;
		%Equivalence of digital image correlation criteria for pattern matching - Bing Pan
		if corrleation_type==3
			%zero mean normalised sum of squared difference
			dG_temp=sum(sum((G_defromed-G_def_mean).^2));
			dG=sqrt(dG_temp);
			for l=1:r
				for j=1:c
					q=q+dfdw{l,j}'*(F(l,j)-Fmean-dF/dG*(G_defromed(l,j)-G_def_mean));
				end
			end
		elseif corrleation_type==2
			%normalised sum of squared difference
			dG_temp=sum(sum((G_defromed-G_def_mean).^2));
			dG=sqrt(dG_temp);
			for l=1:r
				for j=1:c
					q=q+dfdw{l,j}'*(F(l,j)-dF/dG*(G_defromed(l,j)));
				end
			end
		elseif corrleation_type==1
			%zero mean sum of squared difference
			% C_temp=(F-Fmean-(G_defromed-G_def_mean));
			% [r_c,c_c]=size(C_temp);
			% C_temp_vec=reshape(C_temp,[1,r_c*c_c]);

			for l=1:r
				for j=1:c
					q=q+dfdw{l,j}'*(F(l,j)-Fmean-(G_defromed(l,j)-G_def_mean));
				end
			end
		end

		deltaP=-H\q; %solve #4.20
		% determine the new warp function parameters
		% size(deltaP)
		W=updateWarp(P,deltaP);
		% store the warp function parameters
		% P(:,converge+1)=[W(1,3);(W(1,1)-1);W(1,2);W(2,3);W(2,1);(W(2,2)-1)];
		% P(:,converge+1)=[W(1,3);W(2,3)];
		P=Mat2Vec(W);

		% check the stopping criteria value
		check(converge)=sqrt(deltaP(1)^2+(xpmax*deltaP(2))^2+(ypmax*deltaP(3))^2+deltaP(4)^2+(xpmax*deltaP(5))^2+(ypmax*deltaP(6))^2); % #4.24
		% check(converge)=sqrt(deltaP(1)^2+deltaP(2)^2);
		% check(converge)=norm(deltaP);
		% check(converge)=1;

		criteria=0;
		for l=1:r
			for j=1:c
				dx=XX(j)-x0;
				dy=YY(l)-y0;
				% might change
				% xp1(l,j)=x0+dx*(1+P(2,converge))+P(3,converge)*dy+P(1,converge);
				% yp1(l,j)=y0+dy*(1+P(6,converge))+P(5,converge)*dx+P(4,converge);
				% xp1(l,j)=x0+P(1,converge)+dx;
				% yp1(l,j)=y0+P(2,converge)+dy;
				temp=WarpFunc(dx,dy,P);
				xp1(l,j)=x0+temp(1);
				yp1(l,j)=y0+temp(2);
				if corrleation_type==3
					criteria=criteria+((F(l,j)-Fmean)/dF -(G_defromed(l,j)-G_def_mean)/dG)^2;
				elseif corrleation_type==2
					criteria=criteria+((F(l,j))/dF -(G_defromed(l,j))/dG)^2;
				elseif corrleation_type==1
					criteria=criteria+((F(l,j)-Fmean) -(G_defromed(l,j)-G_def_mean))^2;
				end
						
			end
		end

		if (check(converge)<0.0001)&&(criteria<0.001)
			flag=0;
		elseif (converge>2000)
			flag=0;
		end

		% if the subsets are to be considered converged then determine the correlation coefficient #4.10
		if (flag==0)
			% criteria=0;
			% for l=1:r
			% 	for j=1:c
			% 		dx=XX(j)-x0;
			% 		dy=YY(l)-y0;
			% 		% might change
			% 		% xp1(l,j)=x0+dx*(1+P(2,converge))+P(3,converge)*dy+P(1,converge);
			% 		% yp1(l,j)=y0+dy*(1+P(6,converge))+P(5,converge)*dx+P(4,converge);
			% 		% xp1(l,j)=x0+P(1,converge)+dx;
			% 		% yp1(l,j)=y0+P(2,converge)+dy;
			% 		temp=WarpFunc(dx,dy,P);
			% 		xp1(l,j)=x0+temp(1);
			% 		yp1(l,j)=y0+temp(2);
			% 		if corrleation_type==3
			% 			criteria=criteria+((F(l,j)-Fmean)/dF -(G_defromed(l,j)-G_def_mean)/dG)^2;
			% 		elseif corrleation_type==2

			% 		elseif corrleation_type==1
			% 			criteria=criteria+((F(l,j)-Fmean) -(G_defromed(l,j)-G_def_mean))^2;
			% 		end
							
			% 	end
			% end
			Corr_out=criteria;
			% P_final=[P(1,converge+1),P(2,converge+1),P(3,converge+1),P(4,converge+1),P(5,converge+1),P(6,converge+1)];
			P_final=P;
			% fprintf('Number of iterations: %d', converge);
		end
		% P_out(converge,:)=[P(1,converge+1),P(2,converge+1),P(3,converge+1),P(4,converge+1),P(5,converge+1),P(6,converge+1)];
		% P_out(converge,:)=[P(1,converge+1),P(2,converge+1)];
		P_store(:,converge)=P;
		converge=converge+1;
	end
end

%change
function out=Warp(p,d)
	% determine the update to the warp function parameters 
	a=[(1+p(2)), p(3), p(1);
		p(5), (p(6)+1), p(4);
		0 0 1];
	b=[(1+d(2)), d(3), d(1);
		d(5), (d(6)+1), d(4);
		0 0 1];
	out=a*inv(b); % #4.23
end

function out=Warp2(p,d)
	% determine the update to the warp function parameters 
	a=[1, 0, p(1);
		0, 1, p(2);
		0 0 1];
	b=[1, 0, d(1);
		0, 1, d(2);
		0 0 1];
	out=a*inv(b); % #4.23
end

function out=updateWarp(p,d)
	a=WarpMat(p);
	b=WarpMat(d');
	out=a*inv(b); % #4.23
end
%resources
%Advancement of Optical Methods in Experimental Mechanics, Volume 3
%http://www.ncorr.com/index.php/dic-algorithms#3_5
% http://geoserver.ing.puc.cl/info/docencia/ice1603/DIC/Correlation_Tracking_Guide_2010.htm  - DIC tracking program for matlab