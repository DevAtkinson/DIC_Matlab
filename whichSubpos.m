function out=whichSubpos(subpos,stepsize,y,x)
	[r,c]=size(subpos);
	for i=1:r
		for j=1:c
			if (x>=subpos{i,j}.coords(1))&(x<subpos{i,j}.coords(3))&(y>=subpos{i,j}.coords(2))&(y<subpos{i,j}.coords(4))
				out=[i,j];
			end
		end
	end
end