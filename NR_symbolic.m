function NR_symbolic(A,B,save_name,PP)
	%function to create the functions that will give the symbolic jacobian and hessian of the G values
	syms dx dy X Y;
	P = sym('P_%d', size(PP));
		
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

	
	a=reshape(coef,[4,4]);
	% x_dec=mod(xp,1);
	% y_dec=mod(yp,1);
	x_dec=xp-floor(xp);
	y_dec=yp-floor(yp);
	G=[1, x_dec, x_dec^2, x_dec^3]*a*[1; y_dec; y_dec^2; y_dec^3];

	J=jacobian(G,P);
	H=hessian(G,P);
	[folder,name_short,~] = fileparts(save_name);
	name1=fullfile(folder,strcat(name_short,'_jac'));
	name2=fullfile(folder,strcat(name_short,'_hes'));
	name3=fullfile(folder,strcat(name_short,'_pos'));
	warning('off','all');
	matlabFunction(J,'File',name1,'Optimize',true,'Vars',{coef,P,dx,dy,X,Y});
	matlabFunction(H,'File',name2,'Optimize',true,'Vars',{coef,P,dx,dy,X,Y});
	matlabFunction(xp,yp,'File',name3,'Optimize',true,'Vars',{P,dx,dy,X,Y},'Output',{'xp','yp'});
	warning('on','all');
	correctFloor(strcat(name1,'.m'));
	correctFloor(strcat(name2,'.m'));
end