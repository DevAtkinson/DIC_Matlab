function [subpos,mask_subsets,valid_subsets]=mask2subsets(mask,subsize,stepsize)
	[r,c]=size(mask);
	% check=0;
	% check2=0;
	% subcheck=zeros([subsize,subsize]);
	% r_pos=ceil(subsize/2);
	% while check==0
	% 	c_pos=ceil(subsize/2);
	% 	check2==0;
	% 	while check2==0
	% 		for i=1:subsize
	% 			for j=1:subsize
	% 				r_pos+(i-ceil(subsize/2))
	% 				c_pos+(j-ceil(subsize/2))
	% 				subcheck(i,j)=mask(r_pos+(i-ceil(subsize/2)),c_pos+(j-ceil(subsize/2)));
	% 			end
	% 		end

	% 		c_pos=c_pos+1;
	% 		if (sum(sum(subcheck))==subsize*subsize)
	% 			check=1;
	% 			check2=1;
	% 			coords=[r_pos,c_pos-1];
	% 		elseif (c_pos>c)
	% 			check2=1;
	% 		end
	% 	end
	% 	r_pos=r_pos+1;
	% end

	[rmin, cmin, rmax, cmax]=maskOutline(mask)
	r_pos=floor(subsize/2)+rmin;
	c_pos=floor(subsize/2)+cmin;
	coords=[r_pos,c_pos];
	afterr=floor(((rmax-rmin+1)-subsize)/stepsize)+1;
	afterc=floor(((cmax-cmin+1)-subsize)/stepsize)+1;
	for iii=1:stepsize %change position of subsets
		fprintf('%d / %d \n', iii, stepsize);
		for jjj=1:stepsize
			count=0;
			center=[coords(1)+(iii-1),coords(2)+(jjj-1)];
			for ii=1:afterr
				for jj=1:afterc
					cent=[center(1)+stepsize*(ii-1),center(2)+stepsize*(jj-1)];
					mask_check=mask((cent(1)-floor(subsize/2)):(cent(1)+floor(subsize/2)),(cent(2)-floor(subsize/2)):(cent(2)+floor(subsize/2)));
					if prod(prod(mask_check))==1
						count=count+1;
					end
				end
			end
			number(iii,jjj)=count;
		end
	end

	mask_subsets=zeros(size(mask));
	[val,ind]=max(number(:));
	[indr, indc] = ind2sub(size(number),ind);
	center=[coords(1)-indr+1,coords(2)-indc+1];
	valid_subsets=zeros([afterr,afterc]);
	for ii=1:afterr
		for jj=1:afterc
			cent=[center(1)+stepsize*(ii+1),center(2)+stepsize*(jj-1)];
			mask_check=mask((cent(1)-floor(subsize/2)):(cent(1)+floor(subsize/2)),(cent(2)-floor(subsize/2)):(cent(2)+floor(subsize/2)));
			subpos{ii,jj}.coords=[cent(1)-floor(subsize/2),cent(2)-floor(subsize/2),cent(1)+floor(subsize/2),cent(2)+floor(subsize/2)];
			if prod(prod(mask_check))==1
				subpos{ii,jj}.include=1;
				mask_subsets(subpos{ii,jj}.coords(1):subpos{ii,jj}.coords(3),subpos{ii,jj}.coords(2):subpos{ii,jj}.coords(4))=ones(subsize,subsize);
				valid_subsets(ii,jj)=1;
			else
				subpos{ii,jj}.include=0;
			end
		end
	end
			

	% for i=-ceil(subsize/2):1:(subsize-ceil(subsize/2))
	% 	fprintf('%d / %d \n', i, ceil(subsize/2));
	% 	for j=-ceil(subsize/2):1:(subsize-ceil(subsize/2))
	% 		% toc
	% 		% fprintf('%d / %d %d / %d\n', i, ceil(subsize/2), j, ceil(subsize/2));
	% 		% tic
	% 		count=0;
	% 		center=[coords(1)-i,coords(2)-j];
	% 		beforec=floor((center(2)-cmin-floor(subsize))/stepsize);
	% 		afterc=floor((cmax-center(2)-floor(subsize))/stepsize);
	% 		beforer=floor((center(1)-rmin-floor(subsize))/stepsize);
	% 		afterr=floor((rmax-center(1)-floor(subsize))/stepsize);
	% 		for ii=-beforer:1:afterr
	% 			for jj=-beforec:1:afterc
	% 				cent=[center(1)+ii*stepsize,center(2)+jj*stepsize];
	% 				for iii=1:subsize
	% 					for jjj=1:subsize
	% 						subcheck(iii,jjj)=mask(cent(1)+(iii-ceil(subsize/2)),cent(2)+(jjj-ceil(subsize/2)));
	% 					end
	% 				end
	% 				if (sum(sum(subcheck))==subsize*subsize)
	% 					count=count+1;
	% 				end
	% 			end
	% 		end
	% 		number(i+ceil(subsize/2)+1,j+ceil(subsize/2)+1)=count;
	% 	end
	% end











	% r_pos=ceil(subsize/2);
	% c_pos=ceil(subsize/2);
	% coords=[r_pos,c_pos];
	% % tic

	% for i=-ceil(subsize/2):1:(subsize-ceil(subsize/2))
	% 	fprintf('%d / %d \n', i, ceil(subsize/2));
	% 	for j=-ceil(subsize/2):1:(subsize-ceil(subsize/2))
	% 		% toc
	% 		% fprintf('%d / %d %d / %d\n', i, ceil(subsize/2), j, ceil(subsize/2));
	% 		% tic
	% 		count=0;
	% 		center=[coords(1)-i,coords(2)-j];
	% 		beforec=floor((center(2)-1-floor(subsize))/stepsize);
	% 		afterc=floor((c-center(2)-floor(subsize))/stepsize);
	% 		beforer=floor((center(1)-1-floor(subsize))/stepsize);
	% 		afterr=floor((r-center(1)-floor(subsize))/stepsize);
	% 		for ii=-beforer:1:afterr
	% 			for jj=-beforec:1:afterc
	% 				cent=[center(1)+ii*stepsize,center(2)+jj*stepsize];
	% 				for iii=1:subsize
	% 					for jjj=1:subsize
	% 						subcheck(iii,jjj)=mask(cent(1)+(iii-ceil(subsize/2)),cent(2)+(jjj-ceil(subsize/2)));
	% 					end
	% 				end
	% 				if (sum(sum(subcheck))==subsize*subsize)
	% 					count=count+1;
	% 				end
	% 			end
	% 		end
	% 		number(i+ceil(subsize/2)+1,j+ceil(subsize/2)+1)=count;
	% 	end
	% end
	% [val,ind]=max(number(:));
	% [indr, indc] = ind2sub(size(number),ind);
	% indr=indr-ceil(subsize/2)-1;
	% indc=indc-ceil(subsize/2)-1;

	% center=[coords(1)-indr,coords(2)-indc];
	% beforec=floor((center(2)-1-floor(subsize))/stepsize);
	% afterc=floor((c-center(2)-floor(subsize))/stepsize);
	% beforer=floor((center(1)-1-floor(subsize))/stepsize);
	% afterr=floor((r-center(1)-floor(subsize))/stepsize);
	% mask_subsets=zeros(size(mask));
	% for ii=-beforer:1:afterr
	% 	for jj=-beforec:1:afterc
	% 		cent=[center(1)+ii*stepsize,center(2)+jj*stepsize];
	% 		for iii=1:subsize
	% 			for jjj=1:subsize
	% 				subcheck(iii,jjj)=mask(cent(1)+(iii-ceil(subsize/2)),cent(2)+(jjj-ceil(subsize/2)));
	% 			end
	% 		end
	% 		subpos{ii,jj}.coords=[cent(1)-floor(subsize/2),cent(2)-floor(subsize/2),cent(1)+floor(subsize/2),cent(2)+floor(subsize/2)];
	% 		mask_subsets(subpos{ii,jj}.coords(1):subpos{ii,jj}.coords(3),subpos{ii,jj}.coords(2):subpos{ii,jj}.coords(4))=ones(subsize,subsize);
	% 		if (sum(sum(subcheck))==subsize*subsize)
	% 			subpos{ii,jj}.include=1;
	% 		else
	% 			subpos{ii,jj}.include=0;
	% 		end
	% 	end
	% end
end