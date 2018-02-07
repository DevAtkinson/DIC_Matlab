function NR2(subsize,numP,A,X)
	% F=F_in(1:subsize,1:subsize);
	F = sym('F_%d_%d', [subsize subsize]);
	coef = sym('A_%d_%d_%d', [subsize subsize 16]);
	G = sym('G_%d_%d', [subsize subsize]);
	P = sym('P_%d', [numP,1])
	subcord=sym('subcord_%d', [2,1]);
	% syms d	W=A*X;
	% [r_g,c_g]=size(G);
	% [Xmesh,Ymesh]=meshgrid(1:1:subsize,1:1:subsize);
	% Int = (Xmesh',Ymesh',G','cubic')


	% XX=subcord(2):subcord(2)+subsize-1;										% x values over subset
	% YY=subcord(1):subcord(1)+subsize-1;										% y values over subset
	XX=1:subsize;										% x values over subset
	YY=1:subsize;										% y values over subset
	X=repmat(XX,subsize,1);
	Y=repmat(YY',1,subsize);
	% x0=subsize/2+subcord(2);													% x value for subset centre
	% y0=subsize/2+subcord(1);	
	x0=subsize/2;													% x value for subset centre
	y0=subsize/2;	
	% Fmean=mean(mean(F));
	% dF_temp=0;
	dx=X-x0;
	dy=Y-y0;


	% Ftemp=F-ones([subsize, subsize])*Fmean;
	% dF_temp=sum(sum(Ftemp.^2));

	% dF=sqrt(dF_temp);


	
	% dx=reshape(dx,[subsize*subsize,1]);
	% dy=reshape(dy,[subsize*subsize,1]);
	tic
	for i=1:subsize
		for j=1:subsize
			xp(i,j)=P(1)+P(3).*dy(i,j)+dx(i,j).*(P(2)+1.0)+X(i,j);
			yp(i,j)=P(4)+P(5).*dx(i,j)+dy(i,j).*(P(6)+1.0)+Y(i,j);
		end
	end
	toc0=toc
	tic
	% G=zeros([subsize,subsize]);
	for i=1:subsize
		fprintf('interp %d \n', i);
		for j=1:subsize
			a=reshape(coef(i,j,:),[4,4]);
			x_dec=mod(xp(i,j),1);
			y_dec=mod(yp(i,j),1);
			G(i,j)=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];
		end
	end
	toc1=toc
	% tic
	% fprintf('1\n');
	% Gmean=mean(mean(G));
	% fprintf('2\n');
	% Gtemp=G-ones([subsize, subsize])*Gmean;
	% fprintf('3\n');
	% dG_temp=sum(sum(Gtemp.^2));
	% fprintf('4\n');
	% dG=sqrt(dG_temp);
	% toc2=toc
	tic
	% S=0;
	% for i=1:subsize
	% 	fprintf('S %d \n', i);
	% 	for j=1:subsize
	% 		S=S+((F(i,j)-Fmean)/dF -(G(i,j)-Gmean)/dG)^2;
	% 	end
	% end
	S=sum(sum((F-G).^2))/sum(sum(F.^2));
	toc
	% fid = fopen('myfile.txt', 'w');
	% fwrite(fid, char(S), 'char');
	% fclose(fid);
	% matlabFunction(S,'File','SPaper','Optimize',false);
	toc3=toc
	% save('NRS.mat','S');
	tic
	% S=simplify(S);
	% toc
	for i=1:subsize
		for j=1:subsize
			Jtemp{i,j}=jacobian(G(i,j),P);
		end
	end
	toc
	for i=1:subsize
		for j=1:subsize
			for k=1:numP
				Jtemp2{k}(i,j)=Jtemp{i,j}(k);
			end
		end
	end
	toc
	% Jtemp=jacobian(G,P);
	% Jcoef=-2/(sum(sum(F.^2)))*(sum(sum(F-G)));
	Jcoef=-2/(sum(sum(F.^2)));
	for i=1:numP
		J(i)=Jcoef*(sum(sum((F-G)*Jtemp2{i})));
	end

	% J=jacobian(S,P);
	jactoc=toc
	% J(1,1)=diff(S,P(1))
	% J(2,1)=diff(S,P(2))
	% J(3,1)=diff(S,P(2))
	% J(4,1)=diff(S,P(4))
	% J(5,1)=diff(S,P(5))
	% J(6,1)=diff(S,P(6))
	% save('NRfirst.mat','J');
	% fid = fopen('myfile2.txt', 'w');
	% fwrite(fid, char(J), 'char');
	% fclose(fid);
	matlabFunction(J,'File','JacobianPaper','Optimize',false,'Vars',{F,coef,P});
	madeFile=toc

	tic
	Hcoef=2/(sum(sum(F.^2)));
	for i=1:numP
		for j=1:numP
			H(i,j)=Hcoef*sum(sum(Jtemp{i}.*Jtemp{j}));
		end
	end


	% H=hessian(S,P);
	toc
	% for i=1:numP
	% 	for j=1:numP
	% 		H(i,j)=diff(diff(S,P(i)),P(j))
	% 	end
	% end
	% fid = fopen('myfile3.txt', 'w');
	% fwrite(fid, char(H), 'char');
	% fclose(fid);
	matlabFunction(H,'File','HessianPaper','Optimize',false,'Vars',{F,coef,P});
	% save('NRsecond.mat','H');
	toc


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

	% % interpolate the investigated subset to obtain the subset used for comparison purposes
	% G_defromed=InterpFunc(xp,yp);



end