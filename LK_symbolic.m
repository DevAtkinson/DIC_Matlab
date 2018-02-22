function LK_symbolic(A,X,save_name,P_in,distortion_model,vers)
	syms dx dy x0 y0 dxmax dymax ud vd uu vu Xd Yd Xu Yu cx cy;
	k = sym('k_%d', [1 3]);
	% determine warp function equations
	% W=A*X;
	[r,c]=size(A);
	Aold=A;
	AA = sym('A%d%d', [r c]);
	% checking=has(W(1),P1);
	P = sym('P_%d', size(P_in));
	k = sym('k_%d', [1 3]);
		
	% subs(A,P1,P(1))
	for i=1:max(size(P_in))
		A=subs(A,P_in(i),P(i));
	end
	% Wnew=A*X;
	W=A*X;
	Phas=sym('P_%d',[1,(4)*c]);
	count=1;
	for i=1:((2)*c)
		for j=1:((4)*c)
			if (j==i)&sum(has(W(1:2),P(j)))
				dP(:,count)=diff(W,P(j));
				count=count+1;
			end
		end
	end
	% dPcheck


	% count=1;
	% for i=1:((2)*c)
	% 	if (i==1)&sum(has(W(1:2),P1))
	% 		dP(:,count)=diff(W,P1);
	% 	elseif (i==2)&sum(has(W(1:2),P2))
	% 		dP(:,count)=diff(W,P2);
	% 	elseif (i==3)&sum(has(W(1:2),P3))
	% 		dP(:,count)=diff(W,P3);
	% 	elseif (i==4)&sum(has(W(1:2),P4))
	% 		dP(:,count)=diff(W,P4);
	% 	elseif (i==5)&sum(has(W(1:2),P5))
	% 		dP(:,count)=diff(W,P5);
	% 	elseif (i==6)&sum(has(W(1:2),P6))
	% 		dP(:,count)=diff(W,P6);
	% 	elseif (i==7)&sum(has(W(1:2),P7))
	% 		dP(:,count)=diff(W,P7);
	% 	elseif (i==8)&sum(has(W(1:2),P8))
	% 		dP(:,count)=diff(W,P8);
	% 	elseif (i==9)&sum(has(W(1:2),P9))
	% 		dP(:,count)=diff(W,P9);
	% 	elseif (i==10)&sum(has(W(1:2),P10))
	% 		dP(:,count)=diff(W,P10);
	% 	elseif (i==11)&sum(has(W(1:2),P11))
	% 		dP(:,count)=diff(W,P11);
	% 	elseif (i==12)&sum(has(W(1:2),P12))
	% 		dP(:,count)=diff(W,P12);
	% 	elseif (i==13)&sum(has(W(1:2),P13))
	% 		dP(:,count)=diff(W,P13);
	% 	elseif (i==14)&sum(has(W(1:2),P14))
	% 		dP(:,count)=diff(W,P14);
	% 	elseif (i==15)&sum(has(W(1:2),P15))
	% 		dP(:,count)=diff(W,P15);
	% 	elseif (i==16)&sum(has(W(1:2),P16))
	% 		dP(:,count)=diff(W,P16);
	% 	else
	% 		count=count-1;
	% 	end
	% 	count=count+1;
	% end
	% dP

	for i=1:2
		for j=1:(count-1)
			% fprintf('%d %d\n', i,j)
			out(i,j)=dP(i,j);
		end
	end

	for i=1:2
		Wout(i)=W(i);
	end

	% for i=1:2
	% 	for j=1:c
	% 		if has(A(i,j),P1)
				

	% 	end
	% end

	for i=1:2
		for j=1:c
			eqn=A(i,j)==AA(i,j)
			for l=1:max(size(P_in))
				if has(eqn,P(l))
					M(l)=solve(eqn,P(l));
				end
			end
		end
	end

	% for i=1:2
	% 	for j=1:c
	% 		eqn=A(i,j)==AA(i,j)
	% 		if has(eqn,P1)
	% 			M(1)=solve(eqn,P1);
	% 		elseif has(eqn,P2)
	% 			M(2)=solve(eqn,P2);
	% 		elseif has(eqn,P3)
	% 			M(3)=solve(eqn,P3);
	% 		elseif has(eqn,P4)
	% 			M(4)=solve(eqn,P4);
	% 		elseif has(eqn,P5)
	% 			M(5)=solve(eqn,P5);
	% 		elseif has(eqn,P6)
	% 			M(6)=solve(eqn,P6);
	% 		elseif has(eqn,P7)
	% 			M(7)=solve(eqn,P7);
	% 		elseif has(eqn,P8)
	% 			M(8)=solve(eqn,P8);
	% 		elseif has(eqn,P9)
	% 			M(9)=solve(eqn,P9);
	% 		elseif has(eqn,P10)
	% 			M(10)=solve(eqn,P10);
	% 		elseif has(eqn,P11)
	% 			M(11)=solve(eqn,P11);
	% 		elseif has(eqn,P12)
	% 			M(12)=solve(eqn,P12);
	% 		elseif has(eqn,P13)
	% 			M(13)=solve(eqn,P13);
	% 		elseif has(eqn,P14)
	% 			M(14)=solve(eqn,P14);
	% 		elseif has(eqn,P15)
	% 			M(15)=solve(eqn,P15);
	% 		elseif has(eqn,P16)
	% 			M(16)=solve(eqn,P16);
	% 		end
	% 	end
	% end

	%distortion model
	if distortion_model==1
		xu=X*k(1)*((X-cx)^2+(Y-cy)^2);
		yu=Y*k(1)*((X-cx)^2+(Y-cy)^2);
	elseif distortion_model==2
		xu=X*(k(1)*((X-cx)^2+(Y-cy)^2)+k(2)*((X-cx)^2+(Y-cy)^2)^2);
		yu=Y*(k(1)*((X-cx)^2+(Y-cy)^2)+k(2)*((X-cx)^2+(Y-cy)^2)^2);
	elseif distortion_model==3
		xu=X*(k(1)*((X-cx)^2+(Y-cy)^2)+k(2)*((X-cx)^2+(Y-cy)^2)^2+k(3)*((X-cx)^2+(Y-cy)^2)^3);
		yu=Y*(k(1)*((X-cx)^2+(Y-cy)^2)+k(2)*((X-cx)^2+(Y-cy)^2)^2+k(3)*((X-cx)^2+(Y-cy)^2)^3);
	end

	eqxd=subs(xu,X,Xd);
	eqxd=subs(eqxd,Y,Yd);
	eqxu=subs(xu,X,Xu);
	eqxu=subs(eqxu,Y,Yu);

	eqyd=subs(yu,X,Xd);
	eqyd=subs(eqyd,Y,Yd);
	eqyu=subs(yu,X,Xu);
	eqyu=subs(eqyu,Y,Yu);

	uu=P(1)-(eqxd-eqxu);
	vu=P(4)-(eqyd-eqyu);
	uu=subs(uu,Xd,(Xu+P(1)));
	uu=subs(uu,Yd,(Yu+P(4)));

	vu=subs(vu,Xd,(Xu+P(1)));
	vu=subs(vu,Yd,(Yu+P(4)));


	% W_func=matlabFunction(W);
	% out_func=matlabFunction(out);
	% AA=matlabFunction(A);
	% [folder,name_short,~] = fileparts(save_name);
	% name1=fullfile(folder,'WarpFunc'); dx dy P
	% name2=fullfile(folder,'dPFunc'); dx dy -P
	% name3=fullfile(folder,'WarpMat'); P
	% name4=fullfile(folder,'Mat2Vec'); AA

	[folder,name_short,~] = fileparts(save_name);
	name1=fullfile(folder,strcat(name_short,sprintf('_WarpFunc%d',vers)));
	name2=fullfile(folder,strcat(name_short,sprintf('_dPFunc%d',vers)));
	name3=fullfile(folder,strcat(name_short,sprintf('_WarpMat%d',vers)));
	name4=fullfile(folder,strcat(name_short,sprintf('_Mat2Vec%d',vers)));
	name5=fullfile(folder,strcat(name_short,sprintf('_dist%d',vers)));

	matlabFunction(Wout,'File',name1,'Vars',{dx,dy,P});
	matlabFunction(out,'File',name2,'Vars',{dx,dy});
	matlabFunction(A,'File',name3,'Vars',{P});
	matlabFunction(M,'File',name4,'Vars',{[AA]});
	matlabFunction(uu,vu,'File',name5,'Optimize',true,'Vars',{P,Xu,Yu,cx,cy,k},'Output',{'uu','vu'});

	% if (2)*c==6
	% 	matlabFunction(Wout,'File',name1,'Vars',{dx,dy,P});
	% 	matlabFunction(out,'File',name2,'Vars',{dx,dy});
	% 	matlabFunction(A,'File',name3,'Vars',{[P1 P2 P3 P4 P5 P6]});
	% 	matlabFunction(M,'File',name4,'Vars',{[AA]});
	% elseif (2)*c==8
	% 	matlabFunction(Wout,'File',name1,'Vars',{dx,dy,[P1 P2 P3 P4 P5 P6 P7 P8]});
	% 	matlabFunction(out,'File',name2,'Vars',{dx,dy});
	% 	matlabFunction(A,'File',name3,'Vars',{[P1 P2 P3 P4 P5 P6 P7 P8]});
	% 	matlabFunction(M,'File',name4,'Vars',{[AA]});
	% elseif (2)*c==10
	% 	matlabFunction(Wout,'File',name1,'Vars',{dx,dy,[P1 P2 P3 P4 P5 P6 P7 P8 P9 P10]});
	% 	matlabFunction(out,'File',name2,'Vars',{dx,dy});
	% 	matlabFunction(A,'File',name3,'Vars',{[P1 P2 P3 P4 P5 P6 P7 P8 P9 P10]});
	% 	matlabFunction(M,'File',name4,'Vars',{[AA]});
	% end
					
				
end