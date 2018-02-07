function [P_final,Corr_out,iterations]=NRtracking4(varargin)
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
	% size(coef)
	[r_coef,c_coef,~]=size(coef);
	coef2=zeros([r_coef*4,c_coef*4]);
	for i=1:r_coef
		for j=1:c_coef
			% i
			% j
			% squeeze(reshape(coef(i,j,:),[4,4]))
			coef2((i*4-3):(i*4),(j*4-3):(j*4))=reshape(coef(i,j,:),[4,4]);
		end
	end
	% size(coef2)
	% P_temp=P;
	% clear P;
	% P=[P_temp(1), P_temp(4)];
	testing=0;
	dP_flag=0;
	search_flag=0;
	exit_flag=0;
	artificial_dir_count=0;
	fail_count=0;
	normal_fail_count=0;
	stop_check3=0;
	while flag==0
		count=count+1;
		search_length_flag=0;
		% [G,J,H]=getJac1_eff(coef,P,dy,dx,x0,y0,F,2,subpos,stepsize,coef_shift,coef2);
		% meshcompare(F,G)
		if (search_flag==0)%&&(count>2)
			if  count<5 %(count<0)&&(mod(count,10)==0)%(count>1)&((norm(stepChange)<0.000001))&(dP_flag==0)
				if (count==1)||(count==3)
					dP=[1;0;0;0;0;0];
					testing=testing+1;
					dP_flag=1;
				elseif (count==2)|(count==4)
					dP=[0;0;0;1;0;0];
					testing=testing+1;
					dP_flag=1;
				end
				% [G,Funcval]=getJac1_eff(coef,P,dy,dx,x0,y0,F,1,subpos,stepsize,coef_shift);
				% Funcval_start=Funcval
			else
				[G,J,H]=getJac1_eff(coef,P,dy,dx,x0,y0,F,2,subpos,stepsize,coef_shift,coef2,save_name);
				% [G,J,H,Funcval]=getJac5(coef,P,dy,dx,x0,y0,F);
				% [L,D,C]=AdjustCholesky(H,0.0001);
				% H=C*C';
				% dP=abs(H)\J;
				dP=H\J;
				dP_flag=0;
				% [~,pos_def] = chol(H)
			end
			% P_store(count,:)=P;
			dP_store(count,:)=dP';

			if (dP_flag==1)||(normal_fail_count==1)
				[xa,xb]=GoldenBracket(coef,P,dP,dy,dx,x0,y0,F,1,subpos,stepsize,coef_shift,0.02,2,coef2,save_name);
			else
				[xa,xb]=GoldenBracket(coef,P,dP,dy,dx,x0,y0,F,1,subpos,stepsize,coef_shift,0.1,1,coef2,save_name);
			end	

			[G,Funcval,step_change,Pout]=goldenSection(coef,P,dP,dy,dx,x0,y0,F,choice2,subpos,stepsize,coef_shift,save_name,xa,xb);

			if count>2
				valid_check=Funcval/Funcval_store(count-1)<1.2;
				% stop_check2=Funcval_store(count)/Funcval_store(count-1)>0.995;
			else
				% stop_check2=0;
				valid_check=1;
			end
			if valid_check==1%Funcval<Funcval_store(count-1)
				P=Pout;
				P_store(count,:)=P;
				% iterations=i
				stepChange=step_change;
				% norm_step=norm(stepChange);
				stepChange_store(count,:)=stepChange;
				Funcval_store(count)=Funcval;
				artificial_dir_count=0;
				fail_count=0;
				normal_fail_count=0;
				% Funcval
			else
				count=count-1;
				normal_fail_count=normal_fail_count+1;
				if normal_fail_count==1
					search_flag=0;
				elseif normal_fail_count==2
					fail_count=fail_count+1;
					search_flag=1;
				end
				% stepChange=[999,999,999,999,999,999];
				% stepChange_store(count,:)=stepChange;
				% norm_step=norm(stepChange);
				% Funcval_store(count)=Funcval_store(count-1);
				
			end



			% P=P+(dP')*(xa+xb)/2;
			% % iterations=i
			% stepChange=(dP')*(xa+xb)/2
			% norm_step=norm(stepChange);
			% stepChange_store(count,:)=stepChange;

			% [G,Funcval]=getJac1_eff(coef,P,dy,dx,x0,y0,F,1,subpos,stepsize,coef_shift);
			% % [G,J,H,Funcval]=getJac5(coef,P,dy,dx,x0,y0,F);
			% Funcval_store(count)=Funcval;
			% Funcval
		elseif (search_flag==1)%||(count<3)
			dP_temp_store=0;
			% fprintf('Starting direction check\n');
			for j=1:6
				if j==1
					dP=[1; 0; 0; 0; 0; 0];
				elseif j==2
					dP=[0; 0; 0; 1; 0; 0];
				elseif j==3
					dP=[0; 0.05; 0; 0; 0; 0];
				elseif j==4
					dP=[0; 0; 0; 0; 0.05; 0];
				elseif j==5
					dP=[0; 0; 0.05; 0; 0; 0];
				elseif j==6
					dP=[0; 0; 0; 0; 0; 0.05];
				end
				
				[G,Funcval,step_change,Pout]=goldenSection(coef,P,dP,dy,dx,x0,y0,F,choice2,subpos,stepsize,coef_shift,save_name,xa,xb);
				Funcval_check(j)=Funcval;
				stepChange_check(j,:)=step_change;
				P_check(j,:)=Pout;
			end
			% Funcval_check
			[val,minimum]=min(Funcval_check);
			normal_fail_count=0;
			% stepChange_check(minimum,:)
			if count==1
				Funcval_store(count)=val;
				stepChange_store(count,:)=stepChange_check(minimum,:);
				P=P_check(minimum,:);
				P_store(count,:)=P;
				search_flag=0;
			else
				if val<=Funcval_store(count-1);
					Funcval_store(count)=val;
					stepChange_store(count,:)=stepChange_check(minimum,:);
					P=P_check(minimum,:);
					P_store(count,:)=P;
					search_flag=0;
					fail_count=0;
					% artificial_dir_count=artificial_dir_count+1;
				else 
					exit_flag=1;
					count=count-1;
					fail_count=fail_count+1;
					% Funcval_store(count)=Funcval_store(count-1);
					% stepChange_store(count,:)=[999,999,999,999,999,999];
					search_flag=0;
				end
			end

		elseif (search_flag==2)%||(count<3)
			for j=1:2
				if j==1
					dP=[1; 0; 0; 0; 0; 0];
				elseif j==2
					dP=[0; 0; 0; 1; 0; 0];
				end
				
				[G,Funcval,step_change,Pout]=goldenSection(coef,P,dP,dy,dx,x0,y0,F,choice2,subpos,stepsize,coef_shift,save_name,xa,xb);
				P_check(j,:)=Pout;
				stepChange_check(j,:)=step_change;
				Funcval_check(j)=Funcval;
			end
			% Funcval_check
			[val,minimum]=min(Funcval_check);
			normal_fail_count=0;
			% stepChange_check(minimum,:)
			if count==1
				Funcval_store(count)=val;
				stepChange_store(count,:)=stepChange_check(minimum,:);
				P=P_check(minimum,:);
				P_store(count,:)=P;
				search_flag=0;
			else
				if val<=Funcval_store(count-1);
					Funcval_store(count)=val;
					stepChange_store(count,:)=stepChange_check(minimum,:);
					P=P_check(minimum,:);
					P_store(count,:)=P;
					search_flag=0;
					fail_count=0;
					% artificial_dir_count=artificial_dir_count+1;
				else 
					exit_flag=1;
					count=count-1;
					fail_count=fail_count+1;
					% Funcval_store(count)=Funcval_store(count-1);
					% stepChange_store(count,:)=[999,999,999,999,999,999];
					search_flag=0;
				end
			end
		end
		tic
		[G,Funcval]=getJac1_eff(coef,P,dy,dx,x0,y0,F,1,subpos,stepsize,coef_shift,coef2,save_name);
		toc
		% meshcompare(F,G)

		if count>1
			stop_check1=Funcval_store(count)/Funcval_store(count-1)<1.005;
			stop_check2=Funcval_store(count)/Funcval_store(count-1)>0.995;
		else
			stop_check2=0;
			stop_check1=0;
		end
		if (count>90)||((count>10)&&((norm(stepChange_store(count))<1e-7)&&(Funcval_store(count)<0.005)))||((count>20)&&(norm(stepChange_store(count))<1e-10)&&(stop_check1==1)&&(stop_check2==1)&&(stop_check3>5))||(fail_count>20)%||((exit_flag==1)&&count>10)||((artificial_dir_count>10)&&(Funcval_store(count)<0.005)) %||(Funcval<0.001))||(stop_check1==1&&stop_check2==1&&Funcval<0.005&&count>3))%||((norm(stepChange)<0.000001))
			flag=1;
			[val_final,minimum_final]=min(Funcval_store);
			P_final=P_store(minimum_final,:);
			% P_final=P;
			% meshcompare(F,G)
			Corr_out=Funcval_store(minimum_final); 
			% Corr_out=Funcval_store(count); 
			iterations=count;
		% elseif count>10
		% 	if (norm(stepChange_store(count-1,:))<1e-14)&&(norm(stepChange_store(count-2,:))<1e-14)&&(norm(stepChange)<1e-14)%&&(Funcval_store(count)>0.005)
		% 		search_flag=1; %try alternative artificial search directions
		% 	% else
		% 	% 	search_flag=0;
		% 	end
		end
		if ((count>20)&&(norm(stepChange_store(count))<1e-10)&&(stop_check1==1)&&(stop_check2==1))
			search_flag=2;
			stop_check3=stop_check3+1;
		end
	end

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
			[varargout{1},varargoit{2}]=getJac1(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		elseif choice==2
			[varargout{1},varargoit{2},varargout{3}]=getJac1(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		elseif choice==3
			[varargout{1},varargoit{2},varargout{3},varargout{4}]=getJac1(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		end
	elseif choice2==2
		if choice==1
			[varargout{1},varargoit{2}]=getJac2(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		elseif choice==2
			[varargout{1},varargoit{2},varargout{3}]=getJac2(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		elseif choice==3
			[varargout{1},varargoit{2},varargout{3},varargout{4}]=getJac2(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		end
	elseif choice2==3
		if choice==1
			[varargout{1},varargoit{2}]=getJac3(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		elseif choice==2
			[varargout{1},varargoit{2},varargout{3}]=getJac3(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		elseif choice==3
			[varargout{1},varargoit{2},varargout{3},varargout{4}]=getJac3(coef,P,dx,dy,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c);
		end
	end
end

function varargout=getJac1(coef,P,dy,dx,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c)
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

function varargout=getJac2(coef,P,dy,dx,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c)
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

function varargout=getJac3(coef,P,dy,dx,X,Y,F,G,choice,subpos,stepsize,coef_shift,r,c)
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

function [aA,aB]=GoldenBracket(coef,P,dP,dy,dx,X,Y,F,choice,subpos,stepsize,coef_shift,interval,choice2,coef2,save_name)
	% interval=0.05
	if choice2==2
		interval2=interval;
		vector=-1.5:interval2:1.5;
		for i=1:max(size(vector))
			[G,f(i)]=getJac(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
			% if choice==1
			% 	[G,f(i)]=getJac1_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
			% elseif choice==2
			% 	[G,f(i)]=getJac5_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% elseif choice==3
			% 	[G,f(i)]=getJac6_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% end
		end
		% f
		% vector
		% figure
		% plot(vector,f)
		[minval,mini]=min(f);
		if (mini==1)||(mini==max(size(vector)))
			aA=vector(mini);
			aB=vector(mini)+interval2;
		elseif f(mini+1)>f(mini-1)
			aA=vector(mini)-interval2;
			aB=vector(mini);
		else
			aA=vector(mini);
			aB=vector(mini)+interval2;
		end
	elseif choice2==1
		aA=0;
		aB=interval;
	elseif chioce2==3
		interval2=interval;
		vector=-1.5:interval2:1.5;
		for i=1:max(size(vector))
			if choice==1
				[G,f(i),JJ{i},HH{i}]=getJac1_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,3,subpos,stepsize,coef_shift,coef2,save_name);
			elseif choice==2
				[G,f(i)]=getJac5_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			elseif choice==3
				[G,f(i)]=getJac6_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			end
		end
		[minval,mini]=min(f);
		if (mini==1)||(mini==max(size(vector)))
			aA=vector(mini);
			aB=vector(mini)+interval2;
		elseif f(mini+1)>f(mini-1)
			aA=vector(mini)-interval2;
			aB=vector(mini);
		else
			aA=vector(mini);
			aB=vector(mini)+interval2;
		end


		[G,J,H]=getJac1_eff(coef,P,dy,dx,X,Y,F,2,subpos,stepsize,coef_shift,coef2);
		[G,f_0]=getJac1_eff(coef,P,dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2);
		interval2=interval;
		vector=-1.5:interval2:1.5;
		for i=1:max(size(vector))
			if choice==1
				[G,f(i)]=getJac1_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
			elseif choice==2
				[G,f(i)]=getJac5_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			elseif choice==3
				[G,f(i)]=getJac6_eff(coef,(P+vector(i)*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			end
		end
		for i=1:max(size(vector))
			check(i)=(f_0+(J')*(vector(i)*dP'))<f(i);
		end
	end

	% f(mini)=inf;
	% [minval2,mini2]=min(f);
	% aA=min([vector(mini), vector(mini2)]);
	% aB=max([vector(mini), vector(mini2)]);
	% aA=0;
	% aB=interval;
	[G,f1]=getJac(coef,(P+aA*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
	x1=aA;
	[G,f2]=getJac(coef,(P+aB*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
	x2=aB;
	% if choice==1
	% 	[G,f1]=getJac1_eff(coef,(P+aA*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
	% 	% f1=f(xm+aA*U,penalty);
	% 	x1=aA;
	% 	[G,f2]=getJac1_eff(coef,(P+aB*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
	% 	% f2=f(xm+aB*U,penalty);
	% 	x2=aB;
	% elseif choice==2
	% 	[G,f1]=getJac5_eff(coef,(P+aA*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
	% 	% f1=f(xm+aA*U,penalty);
	% 	x1=aA;
	% 	[G,f2]=getJac5_eff(coef,(P+aB*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
	% 	% f2=f(xm+aB*U,penalty);
	% 	x2=aB;
	% elseif choice==3
	% 	[G,f1]=getJac6_eff(coef,(P+aA*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
	% 	% f1=f(xm+aA*U,penalty);
	% 	x1=aA;
	% 	[G,f2]=getJac6_eff(coef,(P+aB*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
	% 	% f2=f(xm+aB*U,penalty);
	% 	x2=aB;
	% end
	d=abs(aB-aA);
	if f1<f2
		x3=aA-d;
	else
		x3=aB+d;
	end
	[G,f3]=getJac(coef,(P+x3*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
	% if choice==1
	% 	[G,f3]=getJac1_eff(coef,(P+x3*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
	% elseif choice==2
	% 	[G,f3]=getJac5_eff(coef,(P+x3*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
	% elseif choice==3
	% 	[G,f3]=getJac6_eff(coef,(P+x3*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
	% end
	% f3=f(xm+x3*U,penalty);
	exit_cond=0;
	count=0;
	d_inc=0;
	if (f3<f1)&&(f3<f2)
		while exit_cond==0 %(f3<min([f1,f2]))
			count=count+1;
			if f1<f2
				x2_temp=x1;
			else
				x2_temp=x2;
			end
			x1_temp=x3;
			x1=min([x1_temp,x2_temp]);
			x2=max([x1_temp,x2_temp]);
			[G,f1]=getJac(coef,(P+x1*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
			[G,f2]=getJac(coef,(P+x2*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
			% if choice==1
			% 	[G,f1]=getJac1_eff(coef,(P+x1*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
			% 	% f1=f(xm+x1*U,penalty);
			% 	[G,f2]=getJac1_eff(coef,(P+x2*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
			% 	% f2=f(xm+x2*U,penalty);
			% elseif choice==2
			% 	[G,f1]=getJac5_eff(coef,(P+x1*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% 	[G,f2]=getJac5_eff(coef,(P+x2*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% elseif choice==3
			% 	[G,f1]=getJac6_eff(coef,(P+x1*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% 	[G,f2]=getJac6_eff(coef,(P+x2*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% end
			% d=abs(x1-x2);
			if count>15
				d=10*d;
				d_inc=d_inc+1;
				count=0;
			end
			if f1<f2
				x3=x1-d;
			else
				x3=x2+d;
			end
			[G,f3]=getJac(coef,(P+x3*dP'),dy,dx,X,Y,F,1,choice,subpos,stepsize,coef_shift,save_name);
			% if choice==1
			% 	[G,f3]=getJac1_eff(coef,(P+x3*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift,coef2,save_name);
			% elseif choice==2
			% 	[G,f3]=getJac5_eff(coef,(P+x3*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% elseif choice==3
			% 	[G,f3]=getJac6_eff(coef,(P+x3*dP'),dy,dx,X,Y,F,1,subpos,stepsize,coef_shift);
			% end
			if f3>min([f1,f2])
				if d_inc==0
					exit_cond=1;
				else
					d=d/10;
					d_inc=d_inc-1;
				end
			end
			% f3=f(xm+x3*U,penalty);
		end
		% until(f3>min([f1,f2]))
	end
	aA=x1;
	aB=x2;
end

function [G,Funcval,step_change,Pout]=goldenSection(coef,P,dP,dy,dx,x0,y0,F,choice2,subpos,stepsize,coef_shift,save_name,xa,xb)
% golden section search
	iterMax=600;
	r=(sqrt(5)-1)/2;
	tol=1e-10;

	Lo=xb-xa;
   	for i=1:iterMax
   		if i==1
   			lambda1=xa+r^2*Lo;
   			lambda2=xa+r*Lo;
   			[G,f1]=getJac(coef,(P+(dP')*lambda1),dy,dx,x0,y0,F,1,choice2,subpos,stepsize,coef_shift,save_name);
   			[G,f2]=getJac(coef,(P+(dP')*lambda2),dy,dx,x0,y0,F,1,choice2,subpos,stepsize,coef_shift,save_name);
			i=i+1; %two function evaluations in this if statement (i starts at one so one is included already)
		end
		if f1>=f2
			xa=lambda1;
			lambda1=lambda2;
			Li=xb-xa;
			lambda2=xa+r*Li;
			evaluat=1;
		elseif f2>f1
			xb=lambda2;
			lambda2=lambda1;
			Li=xb-xa;
			lambda1=xa+r^2*Li;
			evaluat=2;
		else
			disp('There is an error in the golden section code (NaN or inf values)')
		end
		if Li<tol
			break
		else
			if evaluat==1
				f1=f2;
				[G,f2]=getJac(coef,(P+(dP')*lambda2),dy,dx,x0,y0,F,1,choice2,subpos,stepsize,coef_shift,save_name);
			elseif evaluat==2
				f2=f1;
				[G,f1]=getJac(coef,(P+(dP')*lambda1),dy,dx,x0,y0,F,1,choice2,subpos,stepsize,coef_shift,save_name);
			end
		end
	end
	[G,Funcval]=getJac(coef,(P+(dP')*(xa+xb)/2),dy,dx,x0,y0,F,1,choice2,subpos,stepsize,coef_shift,save_name);
	step_change=(dP')*(xa+xb)/2;
	Pout=(P+(dP')*(xa+xb)/2);
end













% function varargout=getJac5_eff(coef,P,dy,dx,X,Y,F,choice,subpos,stepsize,coef_shift)
% 	% zero-mean sum of squared difference
% 	[r,c]=size(F);
% 	xp=zeros([r, c]);
% 	yp=zeros([r, c]);
% 	G=zeros([r, c]);
% 	numP=max(size(P));

% 	% determine the current position of the sample points according to the current estimates of the P parameters
% 	for i=1:r
% 		for j=1:c
% 			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
% 			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
% 		end
% 	end

% 	% determine the G values at the sample points using interpolation
% 	for i=1:r
% 		for j=1:c
% 			a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize-coef_shift(2),:),[4,4]);
% 			x_dec=mod(xp(i,j),1);
% 			y_dec=mod(yp(i,j),1);
% 			if numP==6
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
% 			elseif numP==7
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
% 			elseif numP==8
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
% 			end
% 		end
% 	end
% 	% determine the mean values for the reference and deformed subset
% 	Fmean=mean(mean(F));
% 	Gmean=mean(mean(G));
% 	if choice==1 % if want the correlation coefficient
% 		Corr=sum(sum((F-Fmean - (G-Gmean)).^2));
% 		varargout{1}=G;
% 		varargout{2}=Corr;
% 	elseif choice==2 % if want the Jacobian and Hessian matrices
% 		J=zeros([numP,1]);
% 		H=zeros([numP,numP]);
% 		for i=1:r
% 			for j=1:c
% 				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize,floor(xp(i,j))-subpos.coords(2)+1+stepsize,:),P',dx(i,j),dy(i,j),X,Y));
% 				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize,floor(xp(i,j))-subpos.coords(2)+1+stepsize,:),P',dx(i,j),dy(i,j),X,Y));

% 				J=J+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Jacky');

% 				H=H+(-Jacky')*(-Jacky)+(F(i,j)-Fmean - (G(i,j)-Gmean)).*(-Hess); %if transpose jacky outside of brackets it doesn't work
% 			end
% 		end
% 		J=2*J;
% 		H=2*H;
% 		varargout{1}=G;
% 		varargout{2}=J;
% 		varargout{3}=H;
% 	end
% end

% function varargout=getJac6_eff(coef,P,dy,dx,X,Y,F,choice,subpos,stepsize,coef_shift)
% 	% normalised sum of squared difference
% 	[r,c]=size(F);
% 	xp=zeros([r, c]);
% 	yp=zeros([r, c]);
% 	G=zeros([r, c]);
% 	numP=max(size(P));

% 	% determine the current position of the sample points according to the current estimates of the P parameters
% 	for i=1:r
% 		for j=1:c
% 			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
% 			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
% 		end
% 	end

% 	% determine the G values at the sample points using interpolation
% 	for i=1:r
% 		for j=1:c
% 			a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize-coef_shift(2),:),[4,4]);
% 			x_dec=mod(xp(i,j),1);
% 			y_dec=mod(yp(i,j),1);
% 			if numP==6
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
% 			elseif numP==7
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
% 			elseif numP==8
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
% 			end
% 		end
% 	end
% 	% determine the needed constants for the reference and deformed subset
% 	F2=sqrt(sum(sum(F.^2)));
% 	G2=sqrt(sum(sum(G.^2)));
% 	if choice==1 % if want the correlation coefficient
% 		Corr=sum(sum((F./F2-G./G2).^2));
% 		varargout{1}=G;
% 		varargout{2}=Corr;
% 	elseif choice==2 % if want the Jacobian and Hessian matrices
% 		J=zeros([numP,1]);
% 		H=zeros([numP,numP]);
% 		for i=1:r
% 			for j=1:c
% 				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize,floor(xp(i,j))-subpos.coords(2)+1+stepsize,:),P',dx(i,j),dy(i,j),X,Y));
% 				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize,floor(xp(i,j))-subpos.coords(2)+1+stepsize,:),P',dx(i,j),dy(i,j),X,Y));

% 				J=J+(F(i,j)/F2-G(i,j)/G2).*((-Jacky')/G2);
% 				H=H+(-Jacky'./G2)*(-Jacky./G2) + (F(i,j)/F2-G(i,j)/G2).*(-Hess)/G2; %if transpose jacky outside of brackets it doesn't work
% 			end
% 		end
% 		J=2*J;
% 		H=2*H;
% 		varargout{1}=G;
% 		varargout{2}=J;
% 		varargout{3}=H;
% 	end
% end

% function varargout=getJac1_eff(coef,P,dy,dx,X,Y,F,choice,subpos,stepsize,coef_shift,coef2,save_name)
% 	% zero-mean normalised sum of squared difference
% 	[r,c]=size(F);
% 	xp=zeros([r, c]);
% 	yp=zeros([r, c]);
% 	G=zeros([r, c]);
% 	numP=max(size(P));
% 	% determine the current position of the sample points according to the current estimates of the P parameters
% 	% for i=1:r
% 	% 	for j=1:c
% 	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
% 	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
% 	% 	end
% 	% end
% 	[~,name_short,~] = fileparts(save_name);
% 	funcname=strcat(name_short,'_pos');
% 	funcmeup=str2func(funcname);
% 	[xp,yp]=funcmeup(P',dx,dy,X,Y);
% 	% xp=P(1)+P(3).*dy+dx.*(P(2)+1.0)+X;
% 	% yp=P(4)+P(5).*dx+dy.*(P(6)+1.0)+Y;

% 	% determine the G values at the sample points using interpolation
% 	for i=1:r
% 		for j=1:c
% 			% used=[floor(yp(i,j)),floor(xp(i,j))]
% 			a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2),:),[4,4]);
% 			x_dec=mod(xp(i,j),1);
% 			y_dec=mod(yp(i,j),1);
% 			G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
% 		end
% 	end
	
% 	% for i=1:r
% 	% 	for j=1:c
% 	% 		% used=[floor(yp(i,j)),floor(xp(i,j))]
% 	% 		% a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2),:),[4,4]);
% 	% 		x_dec=mod(xp(i,j),1);
% 	% 		y_dec=mod(yp(i,j),1);
% 	% 		y_ind=floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1);
% 	% 		x_ind=floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2);
% 	% 		X_vec(j,(x_ind*4-3):(x_ind*4))=[1, x_dec, x_dec^2, x_dec^3];
% 	% 		Y_vec((y_ind*4-3):(y_ind*4),i)=[1; y_dec; y_dec^2; y_dec^3];
% 	% 		% if numP==6
% 	% 		% 	G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
% 	% 		% elseif numP==7
% 	% 		% 	G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
% 	% 		% elseif numP==8
% 	% 		% 	G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
% 	% 		% end
% 	% 	end
% 	% end
% 	% % size(X_vec)
% 	% % size(Y_vec)
% 	% G=X_vec*coef2*Y_vec;
% 	% size(G)


% 	% [r_coef,c_coef]=size(coef2);
% 	% X_vec=zeros([r*c,r_coef]);
% 	% Y_vec=zeros([r_coef,c*r]);
% 	% x_dec=mod(xp,1);
% 	% y_dec=mod(yp,1);
% 	% y_ind=floor(yp)-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1);
% 	% x_ind=floor(xp)-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2);
% 	% % y_ind=[2 5; 6 4]
% 	% % r=2
% 	% % c=2
% 	% y_ind1=(reshape(y_ind,1,r*c)-1).*4*r*c+[1:1:r*c];
% 	% y_ind2=y_ind1+1*r*c;
% 	% y_ind3=y_ind1+2*r*c;
% 	% y_ind4=y_ind1+3*r*c;
% 	% % r=2
% 	% % c=2
% 	% % r_coef=160
% 	% % x_ind=[2 3; 4 5]
% 	% x_ind1=((reshape(x_ind,1,r*c)-1).*4+1+([1:1:r*c]-1)*r_coef);
% 	% x_ind2=x_ind1+1;
% 	% x_ind3=x_ind1+2;
% 	% x_ind4=x_ind1+3;

% 	% X_vec(y_ind1)=ones([r,c]);
% 	% X_vec(y_ind2)=x_dec(:);
% 	% X_vec(y_ind3)=x_dec(:).^2;
% 	% X_vec(y_ind4)=x_dec(:).^3;

% 	% Y_vec(x_ind1)=ones([r,c]);
% 	% Y_vec(x_ind2)=y_dec(:);
% 	% Y_vec(x_ind3)=y_dec(:).^2;
% 	% Y_vec(x_ind4)=y_dec(:).^3;
% 	% % size(Y_vec)
% 	% % Y_vec
% 	% GG=X_vec*coef2*Y_vec;
% 	% G=reshape(diag(GG),r,c);



% 	% determine the needed constants for the reference and deformed subset
% 	Fmean=mean(mean(F));
% 	Gmean=mean(mean(G));
% 	F2=sqrt(sum(sum((F-Fmean).^2)));
% 	G2=sqrt(sum(sum((G-Gmean).^2)));
% 	if choice==1 % if want the correlation coefficient
% 		% Corr=0;
% 		% for i=1:r
% 		% 	for j=1:c
% 		% 		Corr=Corr+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).^2;
% 		% 	end
% 		% end
% 		Corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
% 		varargout{1}=G;
% 		varargout{2}=Corr;
% 	elseif choice==2 % if want the Jacobian and Hessian matrices
% 		J=zeros([numP,1]);
% 		H=zeros([numP,numP]);
% 		for i=1:r
% 			for j=1:c
% 				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
% 				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

% 				J=J+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*((-Jacky')./G2);
% 				H=H+(-Jacky'./G2)*(-Jacky./G2) + ((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*(-Hess)./G2; %if transpose jacky outside of brackets it doesn't work
% 			end
% 		end
% 		J=2*J;
% 		H=2*H;
% 		varargout{1}=G;
% 		varargout{2}=J;
% 		varargout{3}=H;
% 	elseif choice==3
% 		Corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
% 		J=zeros([numP,1]);
% 		H=zeros([numP,numP]);
% 		for i=1:r
% 			for j=1:c
% 				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
% 				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

% 				J=J+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*((-Jacky')./G2);
% 				H=H+(-Jacky'./G2)*(-Jacky./G2) + ((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*(-Hess)./G2; %if transpose jacky outside of brackets it doesn't work
% 			end
% 		end
% 		J=2*J;
% 		H=2*H;
% 		varargout{1}=G;
% 		varargout{2}=Corr;
% 		varargout{3}=J;
% 		varargout{4}=H;
% 	end
% end












% function varargout=getJac1_eff(coef,P,dy,dx,X,Y,F,choice,subpos,stepsize,coef_shift,coef2)
% 	% zero-mean normalised sum of squared difference
% 	[r,c]=size(F);
% 	xp=zeros([r, c]);
% 	yp=zeros([r, c]);
% 	G=zeros([r, c]);
% 	numP=max(size(P));
% 	% determine the current position of the sample points according to the current estimates of the P parameters
% 	% for i=1:r
% 	% 	for j=1:c
% 	% 		xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X;
% 	% 		yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y;
% 	% 	end
% 	% end
% 	xp=P(1)+P(3).*dy+dx.*(P(2)+1.0)+X;
% 	yp=P(4)+P(5).*dx+dy.*(P(6)+1.0)+Y;

% 	% determine the G values at the sample points using interpolation
% 	for i=1:r
% 		for j=1:c
% 			% used=[floor(yp(i,j)),floor(xp(i,j))]
% 			a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2),:),[4,4]);
% 			x_dec=mod(xp(i,j),1);
% 			y_dec=mod(yp(i,j),1);
% 			if numP==6
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
% 			elseif numP==7
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
% 			elseif numP==8
% 				G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
% 			end
% 		end
% 	end
	
% 	% for i=1:r
% 	% 	for j=1:c
% 	% 		% used=[floor(yp(i,j)),floor(xp(i,j))]
% 	% 		% a=reshape(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2),:),[4,4]);
% 	% 		x_dec=mod(xp(i,j),1);
% 	% 		y_dec=mod(yp(i,j),1);
% 	% 		y_ind=floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1);
% 	% 		x_ind=floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2);
% 	% 		X_vec(j,(x_ind*4-3):(x_ind*4))=[1, x_dec, x_dec^2, x_dec^3];
% 	% 		Y_vec((y_ind*4-3):(y_ind*4),i)=[1; y_dec; y_dec^2; y_dec^3];
% 	% 		% if numP==6
% 	% 		% 	G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
% 	% 		% elseif numP==7
% 	% 		% 	G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3]+P(7);
% 	% 		% elseif numP==8
% 	% 		% 	G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3].*P(8)+P(7);
% 	% 		% end
% 	% 	end
% 	% end
% 	% % size(X_vec)
% 	% % size(Y_vec)
% 	% G=X_vec*coef2*Y_vec;
% 	% size(G)


% 	% [r_coef,c_coef]=size(coef2);
% 	% X_vec=zeros([r*c,r_coef]);
% 	% Y_vec=zeros([r_coef,c*r]);
% 	% x_dec=mod(xp,1);
% 	% y_dec=mod(yp,1);
% 	% y_ind=floor(yp)-subpos.coords(1)+1+stepsize*coef_shift(3)-coef_shift(1);
% 	% x_ind=floor(xp)-subpos.coords(2)+1+stepsize*coef_shift(3)-coef_shift(2);
% 	% % y_ind=[2 5; 6 4]
% 	% % r=2
% 	% % c=2
% 	% y_ind1=(reshape(y_ind,1,r*c)-1).*4*r*c+[1:1:r*c];
% 	% y_ind2=y_ind1+1*r*c;
% 	% y_ind3=y_ind1+2*r*c;
% 	% y_ind4=y_ind1+3*r*c;
% 	% % r=2
% 	% % c=2
% 	% % r_coef=160
% 	% % x_ind=[2 3; 4 5]
% 	% x_ind1=((reshape(x_ind,1,r*c)-1).*4+1+([1:1:r*c]-1)*r_coef);
% 	% x_ind2=x_ind1+1;
% 	% x_ind3=x_ind1+2;
% 	% x_ind4=x_ind1+3;

% 	% X_vec(y_ind1)=ones([r,c]);
% 	% X_vec(y_ind2)=x_dec(:);
% 	% X_vec(y_ind3)=x_dec(:).^2;
% 	% X_vec(y_ind4)=x_dec(:).^3;

% 	% Y_vec(x_ind1)=ones([r,c]);
% 	% Y_vec(x_ind2)=y_dec(:);
% 	% Y_vec(x_ind3)=y_dec(:).^2;
% 	% Y_vec(x_ind4)=y_dec(:).^3;
% 	% % size(Y_vec)
% 	% % Y_vec
% 	% GG=X_vec*coef2*Y_vec;
% 	% G=reshape(diag(GG),r,c);



% 	% determine the needed constants for the reference and deformed subset
% 	Fmean=mean(mean(F));
% 	Gmean=mean(mean(G));
% 	F2=sqrt(sum(sum((F-Fmean).^2)));
% 	G2=sqrt(sum(sum((G-Gmean).^2)));
% 	if choice==1 % if want the correlation coefficient
% 		% Corr=0;
% 		% for i=1:r
% 		% 	for j=1:c
% 		% 		Corr=Corr+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).^2;
% 		% 	end
% 		% end
% 		Corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
% 		varargout{1}=G;
% 		varargout{2}=Corr;
% 	elseif choice==2 % if want the Jacobian and Hessian matrices
% 		J=zeros([numP,1]);
% 		H=zeros([numP,numP]);
% 		for i=1:r
% 			for j=1:c
% 				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
% 				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

% 				J=J+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*((-Jacky')./G2);
% 				H=H+(-Jacky'./G2)*(-Jacky./G2) + ((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*(-Hess)./G2; %if transpose jacky outside of brackets it doesn't work
% 			end
% 		end
% 		J=2*J;
% 		H=2*H;
% 		varargout{1}=G;
% 		varargout{2}=J;
% 		varargout{3}=H;
% 	elseif choice==3
% 		Corr=sum(sum(((F-Fmean)./F2-(G-Gmean)./G2).^2));
% 		J=zeros([numP,1]);
% 		H=zeros([numP,numP]);
% 		for i=1:r
% 			for j=1:c
% 				Jacky=(JacobianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));
% 				Hess=(HessianStandard(coef(floor(yp(i,j))-subpos.coords(1)+1+stepsize*coef_shift(3),floor(xp(i,j))-subpos.coords(2)+1+stepsize*coef_shift(3),:),P',dx(i,j),dy(i,j),X,Y));

% 				J=J+((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*((-Jacky')./G2);
% 				H=H+(-Jacky'./G2)*(-Jacky./G2) + ((F(i,j)-Fmean)/F2-(G(i,j)-Gmean)/G2).*(-Hess)./G2; %if transpose jacky outside of brackets it doesn't work
% 			end
% 		end
% 		J=2*J;
% 		H=2*H;
% 		varargout{1}=G;
% 		varargout{2}=Corr;
% 		varargout{3}=J;
% 		varargout{4}=H;
% 	end
% end