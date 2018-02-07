function coef=getBicubicValues(F)
	[r,c]=size(F);
	coef=zeros([r-2,c-2,16]);
	for i=2:r-2
		for j=2:c-2
			% coef{i,j}=bicubicCoefficients(F(i-1:i+2,j-1:j+2));
			coef(i,j,:)=reshape(bicubicCoefficients(F(i-1:i+2,j-1:j+2)),[1,1,16]);
		end
	end

end