function [rmin, cmin, rmax, cmax]=maskOutline(mask)
	[r,c]=size(mask);
	for i=1:r
		if sum(mask(i,:))==0
			cmin_temp(i)=NaN;
		else
			cmin_temp(i)=min(find(mask(i,:)));
		end

		if sum(mask(i,:))==0
			cmax_temp(i)=NaN;
		else
			cmax_temp(i)=max(find(mask(i,:)));
		end
	end

	for i=1:c
		if sum(mask(:,i))==0
			rmin_temp(i)=NaN;
		else
			rmin_temp(i)=min(find(mask(:,i)));
		end

		if sum(mask(:,i))==0
			rmax_temp(i)=NaN;
		else
			rmax_temp(i)=max(find(mask(:,i)));
		end
	end
	rmin=min(rmin_temp);
	rmax=max(rmax_temp);
	cmin=min(cmin_temp);
	cmax=max(cmax_temp);

end