function NR_symbolic(A,B,save_name,PP,distortion_model,vers)
	%function to create the functions that will give the symbolic jacobian and hessian of the G values
	syms dx dy X Y k1 k2 k3 ud vd uu vu Xd Yd Xu Yu cx cy;
	P = sym('P_%d', size(PP));
	k = sym('k_%d', [1 3]);
		
	% subs(A,P1,P(1))
	for i=1:max(size(PP))
		A=subs(A,PP(i),P(i));
	end

	XX=A*B;
	xp=XX(1)+X;
	yp=XX(2)+Y;

	coef = sym('A_%d', [1 16]);

	% tic
	% syms dx dy X Y;
	% xp=P(1)+P(3).*dy+dx.*(P(2)+1.0)+X;
	% yp=P(4)+P(5).*dx+dy.*(P(6)+1.0)+Y;

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

	%determine G
	a=reshape(coef,[4,4]);
	% x_dec=mod(xp,1);
	% y_dec=mod(yp,1);
	x_dec=xp-floor(xp);
	y_dec=yp-floor(yp);
	G=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];

	J=jacobian(G,P);
	H=hessian(G,P);
	[folder,name_short,~] = fileparts(save_name);
	name1=fullfile(folder,strcat(name_short,sprintf('_jac%d',vers)));
	name2=fullfile(folder,strcat(name_short,sprintf('_hes%d',vers)));
	name3=fullfile(folder,strcat(name_short,sprintf('_pos%d',vers)));
	name4=fullfile(folder,strcat(name_short,sprintf('_dist%d',vers)));
	warning('off','all');
	matlabFunction(J,'File',name1,'Optimize',true,'Vars',{coef,P,dx,dy,X,Y});
	matlabFunction(H,'File',name2,'Optimize',true,'Vars',{coef,P,dx,dy,X,Y});
	matlabFunction(xp,yp,'File',name3,'Optimize',true,'Vars',{P,dx,dy,X,Y},'Output',{'xp','yp'});
	matlabFunction(uu,vu,'File',name4,'Optimize',true,'Vars',{P,Xu,Yu,cx,cy,k},'Output',{'uu','vu'});
	warning('on','all');
	correctFloor(strcat(name1,'.m'));
	correctFloor(strcat(name2,'.m'));
end