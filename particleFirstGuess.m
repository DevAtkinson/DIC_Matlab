function [Pout,Funcval]=particleFirstGuess(varargin)
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
		case 'stepsize'
			stepsize=varargin{i*2};
		case 'coef_shift'
			coef_shift=varargin{i*2};
		case 'save_name'
			save_name=varargin{i*2};
		case 'algorithm'
			choice2=varargin{i*2};
		case 'lb'
			lb=varargin{i*2};
		case 'ub'
			ub=varargin{i*2};
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


	% options = optimoptions('particle','Display','iter');
	% lb=P-0.2*abs(P);
	% ub=P+0.2*abs(P);
	% tic
	Pout = particleswarm(@(P) getJac(coef,P,dy,dx,x0,y0,F,4,choice2,subpos,stepsize,coef_shift,save_name),6,lb,ub);
	% toc
	[G,Funcval]=getJac(coef,Pout,dy,dx,x0,y0,F,1,choice2,subpos,stepsize,coef_shift,save_name);
	% Funcval
	% Funcval_store(count)=Funcval;
	% stepChange_store(count,:)=Pout-P;
	% P=Pout;
	% P_store(count,:)=P;

end



function varargout=getJac(coef,P,dy,dx,X,Y,F,choice,choice2,subpos,stepsize,coef_shift,save_name)
	% zero-mean normalised sum of squared difference
	[r,c]=size(F);
	xp=zeros([r, c]);
	yp=zeros([r, c]);
	G=zeros([r, c]);
	% numP=max(size(P));

	% determine the current position of the sample points according to the current estimates of the P parameters
	[~,name_short,~] = fileparts(save_name);
	funcname=strcat(name_short,'_pos');
	funcmeup=str2func(funcname);
	[xp,yp]=funcmeup(P',dx,dy,X,Y);

	% determine the G values at the sample points using interpolation
	for i=1:r
		for j=1:c
			% used=[floor(yp(i,j)),floor(xp(i,j))]
			a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2),:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
		end
	end

	% determine the correlation coefficient or jacobian and hessian according to the chosen correlation coefficient equation
	if choice2==1
		if choice==1
			[varargout{1},varargout{2}]=getJac1(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==2
			[varargout{1},varargout{2},varargout{3}]=getJac1(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==3
			[varargout{1},varargout{2},varargout{3},varargout{4}]=getJac1(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==4
			[~,varargout{1}]=getJac1(coef,P,dx,dy,X,Y,F,G,1,subpos,stepsize,coef_shift,r,c,xp,yp);
		end
	elseif choice2==2
		if choice==1
			[varargout{1},varargout{2}]=getJac2(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==2
			[varargout{1},varargout{2},varargout{3}]=getJac2(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==3
			[varargout{1},varargout{2},varargout{3},varargout{4}]=getJac2(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==4
			[~,varargout{1}]=getJac2(coef,P,dx,dy,X,Y,F,G,1,subpos,stepsize,coef_shift,r,c,xp,yp);
		end
	elseif choice2==3
		if choice==1
			[varargout{1},varargout{2}]=getJac3(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==2
			[varargout{1},varargout{2},varargout{3}]=getJac3(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==3
			[varargout{1},varargout{2},varargout{3},varargout{4}]=getJac3(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp);
		elseif choice==4
			[~,varargout{1}]=getJac3(coef,P,dx,dy,X,Y,F,G,1,subpos,stepsize,coef_shift,r,c,xp,yp);
		end
	end
end

function varargout=getJac1(coef,P,dy,dx,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp)
	% zero-mean normalised sum of squared difference

	% determine the needed constants for the reference and deformed subset
	numP=max(size(P));
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));
	F2=sqrt(sum(sum((F-Fmean).^2)));
	G2=sqrt(sum(sum((G-Gmean).^2)));
	if choice==1 % if want the correlation coefficient
		Corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
		varargout{1}=G;
		varargout{2}=Corr;
	elseif choice==2 % if want the Jacobian and Hessian matrices
		J=zeros([numP,1]);
		H=zeros([numP,numP]);
		for i=1:r
			for j=1:c
				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

				J=J+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*((-Jacky')./G2);
				H=H+(-Jacky'./G2)*(-Jacky./G2) + ((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*(-Hess)./G2; %if transpose jacky outside of brackets it doesn't work
			end
		end
		J=2*J;
		H=2*H;
		varargout{1}=G;
		varargout{2}=J;
		varargout{3}=H;
	elseif choice==3
		Corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
		J=zeros([numP,1]);
		H=zeros([numP,numP]);
		for i=1:r
			for j=1:c
				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

				J=J+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*((-Jacky')./G2);
				H=H+(-Jacky'./G2)*(-Jacky./G2) + ((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*(-Hess)./G2; %if transpose jacky outside of brackets it doesn't work
			end
		end
		J=2*J;
		H=2*H;
		varargout{1}=G;
		varargout{2}=Corr;
		varargout{3}=J;
		varargout{4}=H;
	end
end

function varargout=getJac2(coef,P,dy,dx,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp)
	% normalised sum of squared difference
	
	% determine the needed constants for the reference and deformed subset
	numP=max(size(P));
	F2=sqrt(sum(sum(F.^2)));
	G2=sqrt(sum(sum(G.^2)));
	if choice==1 % if want the correlation coefficient
		Corr=sum(sum((F./F2-G./G2).^2));
		varargout{1}=G;
		varargout{2}=Corr;
	elseif choice==2 % if want the Jacobian and Hessian matrices
		J=zeros([numP,1]);
		H=zeros([numP,numP]);
		for i=1:r
			for j=1:c
				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

				J=J+(F(i,j)/F2-G(i,j)/G2).*((-Jacky')/G2);
				H=H+(-Jacky'./G2)*(-Jacky./G2) + (F(i,j)/F2-G(i,j)/G2).*(-Hess)/G2; %if transpose jacky outside of brackets it doesn't work
			end
		end
		J=2*J;
		H=2*H;
		varargout{1}=G;
		varargout{2}=J;
		varargout{3}=H;
	elseif choice==3 % if want the Jacobian and Hessian matrices
		J=zeros([numP,1]);
		H=zeros([numP,numP]);
		for i=1:r
			for j=1:c
				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

				J=J+(F(i,j)/F2-G(i,j)/G2).*((-Jacky')/G2);
				H=H+(-Jacky'./G2)*(-Jacky./G2) + (F(i,j)/F2-G(i,j)/G2).*(-Hess)/G2; %if transpose jacky outside of brackets it doesn't work
			end
		end
		J=2*J;
		H=2*H;
		Corr=sum(sum((F./F2-G./G2).^2));
		varargout{1}=G;
		varargout{2}=Corr;
		varargout{3}=J;
		varargout{4}=H;
	end
end

function varargout=getJac3(coef,P,dy,dx,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c,xp,yp)
	% zero-mean sum of squared difference
	numP=max(size(P));
	% determine the mean values for the reference and deformed subset
	Fmean=mean(mean(F));
	Gmean=mean(mean(G));
	if choice==1 % if want the correlation coefficient
		Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
		varargout{1}=G;
		varargout{2}=Corr;
	elseif choice==2 % if want the Jacobian and Hessian matrices
		J=zeros([numP,1]);
		H=zeros([numP,numP]);
		for i=1:r
			for j=1:c
				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

				J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');

				H=H+(-Jacky')*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess); %if transpose jacky outside of brackets it doesn't work
			end
		end
		J=2*J;
		H=2*H;
		varargout{1}=G;
		varargout{2}=J;
		varargout{3}=H;
	elseif choice==3 % if want the Jacobian and Hessian matrices
		Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
		J=zeros([numP,1]);
		H=zeros([numP,numP]);
		for i=1:r
			for j=1:c
				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

				J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');

				H=H+(-Jacky')*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess); %if transpose jacky outside of brackets it doesn't work
			end
		end
		J=2*J;
		H=2*H;
		varargout{1}=G;
		varargout{2}=Corr;
		varargout{3}=J;
		varargout{4}=H;
	end
end