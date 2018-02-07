function out=Cholesky(L,D)
	A=L*D*L';
	% out=chol(A,'lower');

	[n nn]=size(A);
	for k=1:n
	    A(k,k)=sqrt(A(k,k));
	    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
	    for j=k+1:n
	        A(j:n,j)=A(j:n,j)-A(j,k)*A(j:n,k);
	    end
	end
	out=tril(A);

end

% function L=cholesky (a)
% % Descomposicin Cholesky
% % a: Matriz cuadrada
% [n]=size(a,1);
% [L]=zeros(n,n);
% for k = 1 : n
%    for i = 1 : k-1
%      soma = 0;
%      for j = 1 : i-1
%        soma = soma + a(i,j) * a(k,j);
%      end
%      a(k, i) = (a(k, i) - soma)/a(i, i);
%      if a(k,i)>bound
%      	a(i,i)=a(i,i)*a(k,i)/bound;
%      end
%    end
%    soma = 0;
%    for j = 1 : k -1
%      soma = soma + a(k,j)^2;
%    end
% 	a(k,k) = (a(k,k) - soma)^.5;
% end
% L=tril(a);

% function L=cholesky (a)
% % Descomposicin Cholesky
% % a: Matriz cuadrada
% [n]=size(a,1);
% [L]=zeros(n,n);
% for k = 1 : n
%    for i = 1 : k-1
%      soma = 0;
%      for j = 1 : i-1
%        soma = soma + a(i,j) * a(k,j);
%      end
%      a(k, i) = (a(k, i) - soma)/a(i, i);
%    end
%    soma = 0;
%    for j = 1 : k -1
%      soma = soma + a(k,j)^2;
%    end
% a(k,k) = (a(k,k) - soma)^.5;
% end
% L=tril(a);