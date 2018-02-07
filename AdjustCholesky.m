function [L,D,C]=AdjustCholesky(A,mins)
	[L,D,P]=ldl(A);
	for i=1:max(size(A))
		if abs(D(i,i))<mins
			D(i,i)=mins;
		else
			D(i,i)=abs(D(i,i));
		end
	end
	count=1;
	for i=-(max(size(A))):1:max(size(A))
		modulus(count)=norm(diag(A,i));
		if i~=0
			modulus(count)=modulus(count)/(max(size(A)));
		end
		count=count+1;
	end
	modulus(count)=mins;
	% modulus
	bound=sqrt(max(modulus))

	C=Cholesky(L,D)
	[r,c]=size(C);
	while sum(sum(abs(C)>bound))>0
		check=abs(C)>bound;
		flag=0;
		% for j=1:c
		% 	for i=1:r
		% 		if (check(i,j)==1)&(flag==0)
		% 			fprintf('row: %d col: %d\n',i,j);
		% 			if i==j
		% 				summing=0;
		% 				if j>1
		% 					for jj=1:j-1
		% 						summing=summing+C(j,jj)^2;
		% 					end
		% 				end
		% 				A(j,j)=bound^2+summing;
		% 			else
		% 				summing=0;
		% 				if j>1
		% 					for jj=1:j-1
		% 						summing=summing+C(j,jj)*C(i,jj);
		% 					end
		% 				end
		% 				% C_diagonal_to_change=summing+bound*C(j,j);
		% 				A
		% 				A(j,i)=summing+bound*C(j,j)
		% 				% summing=0;
		% 				% if j>1
		% 				% 	for jj=1:j-1
		% 				% 		summing=summing+C(j,jj)^2;
		% 				% 	end
		% 				% end
		% 				% C_diagonal_to_change=C_diagonal_to_change_temp^2+summing;
		% 			end
		% 			% D(j,j)=C_diagonal_to_change;
		% 			% D(j,j)=D(j,j)*(bound/C(i,j))^2;
		% 			flag=1;

		% 			[L,D,P]=ldl(A);
		% 			for i=1:max(size(A))
		% 				if abs(D(i,i))<mins
		% 					D(i,i)=mins;
		% 				else
		% 					D(i,i)=abs(D(i,i));
		% 				end
		% 			end
		% 			break
		% 		end
		% 	end
		% 	if flag==1
		% 		break
		% 	end
		% end



		for j=1:c
			for i=1:r
				if check(i,j)==1
					fprintf('row: %d col: %d\n',i,j);
					if i==j
						summing=0;
						if j>1
							for jj=1:j-1
								summing=summing+C(j,jj)^2;
							end
						end
						D(j,j)=bound^2+summing;
					else
						summing=0;
						if j>1
							for jj=1:j-1
								summing=summing+C(j,jj)*C(i,jj);
							end
						end
						C_diagonal_to_change_temp=(A(j,i)-summing)/bound;
						summing=0;
						if j>1
							for jj=1:j-1
								summing=summing+C(j,jj)^2;
							end
						end
						C_diagonal_to_change=C_diagonal_to_change_temp^2+summing;
					end
					D(j,j)=C_diagonal_to_change;
					% D(j,j)=D(j,j)*(bound/C(i,j))^2;
					flag=1;
					break
				end
			end
			if flag==1
				break
			end
		end


		% for j=1:c
		% 	for i=1:r
		% 		if check(i,j)==1
		% 			fprintf('row: %d col: %d\n',i,j);
		% 			if i==j
		% 				summing=0;
		% 				if j>1
		% 					for jj=1:j-1
		% 						summing=summing+C(i,jj)^2;
		% 					end
		% 				end
		% 				D(j,j)=bound^2+summing;
		% 			else
		% 				if j>1
		% 					for jj=1:j-1
		% 						summing=summing+C(i,jj)^2;
		% 					end
		% 				end
		% 				C_diagonal_to_change=bound^2+summing;
		% 			end
		% 			D(j,j)=D(j,j)*(bound/C(i,j))^2;
		% 			flag=1;
		% 			break
		% 		end
		% 	end
		% 	if flag==1
		% 		break
		% 	end
		% end
		C=Cholesky(L,D)
	end





	% while max(max(abs(C)))>bound
	% 	% fprintf('Cholesky\n');
	% 	[maximum,ind]=max(abs(C(:)));
	% 	[ind_r,ind_c]=ind2sub(size(C),ind)
	% 	if ind_r==ind_c
	% 		D(ind_c,ind_c)=D(ind_c,ind_c)*bound/C(ind_r,ind_c);
	% 	else
	% 		D(ind_c,ind_c)=D(ind_c,ind_c)*bound/C(ind_r,ind_c);
	% 		% D(ind_c,ind_c)=D(ind_c,ind_c)*C(ind_r,ind_c)/bound;
	% 	end
	% 	% D(ind_c,ind_c)=D(ind_c,ind_c)*1.1;
	% 	C=Cholesky(L,D)

	% end
end




%
% [L,D]=ldlt(A)
%
% This function computes the square root free Cholesky factorization
%
%    A=L*D*L'
%
% where L is a lower triangular matrix with ones on the diagonal, and D
% is a diagonal matrix.  
%
% It is assumed that A is symmetric and postive definite.
%
% Reference: Golub and Van Loan, "Matrix Computations", second edition, 
%            p 137.  
% Author: Brian Borchers (borchers@nmt.edu)
%
function [L,D]=ldlt(A)
	%
	%  Figure out the size of A.
	%
	n=size(A,1);
	%
	% The main loop.  See Golub and Van Loan for details.  
	%
	L=zeros(n,n);
	for j=1:n,
	  if (j > 1),
	    v(1:j-1)=L(j,1:j-1).*d(1:j-1);
	    v(j)=A(j,j)-L(j,1:j-1)*v(1:j-1)';
	    d(j)=v(j);
	    if (j < n),
	      L(j+1:n,j)=(A(j+1:n,j)-L(j+1:n,1:j-1)*v(1:j-1)')/v(j);
	    end;
	  else
	    v(1)=A(1,1);
	    d(1)=v(1);
	    L(2:n,1)=A(2:n,1)/v(1);    
	  end;
	end;
	%
	%  Put d into a matrix.
	%
	D=diag(d);
	%
	%  Put ones on the diagonal of L.
	%
	L=L+eye(n);
end