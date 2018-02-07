function [process_sorted,getindex]=correlationOrder(subpos,selected)
	[r,c]=size(subpos)
	% count=zeros(max([r,c]));
	count=1;
	for i=1:r
		for j=1:c
			dis=max([abs(selected(1)-i), abs(selected(2)-j)]);
			if subpos{i,j}.include==1
				A(i,j)=dis;
			else
				A(i,j)=inf;
			end
			count=count+1;
		end
	end

	count=1;
	for i=1:r
		for j=1:c
			dis=max([abs(selected(1)-i), abs(selected(2)-j)]);
			if subpos{i,j}.include==1
				if (i~=1)&(A(i-1,j)<A(i,j))
					process(count,:)=[dis,i,j,(i-1),j];
				elseif (j~=1)&(A(i,j-1)<A(i,j))
					process(count,:)=[dis,i,j,i,(j-1)];
				elseif (i~=r)&(A(i+1,j)<A(i,j))
					process(count,:)=[dis,i,j,(i+1),j];
				elseif (j~=c)&(A(i,j+1)<A(i,j))
					process(count,:)=[dis,i,j,i,(j+1)];
				elseif (j~=c)&(i~=r)&(A(i+1,j+1)<A(i,j))
					process(count,:)=[dis,i,j,(i+1),(j+1)];
				elseif (j~=c)&(i~=1)&(A(i-1,j+1)<A(i,j))
					process(count,:)=[dis,i,j,(i-1),(j+1)];
				elseif (j~=1)&(i~=r)&(A(i+1,j-1)<A(i,j))
					process(count,:)=[dis,i,j,(i+1),(j-1)];
				elseif (j~=1)&(i~=1)&(A(i-1,j-1)<A(i,j))
					process(count,:)=[dis,i,j,(i-1),(j-1)];
				end
				count=count+1;
			end
		end
	end



	process_sorted=sortrows(process);
	process_sorted(1,:)=[0,selected(1),selected(2),0,0];

	for i=1:max(size(process_sorted))
		count=1;
		for j=1:max(size(process_sorted))
			if (process_sorted(i,2)==process_sorted(j,4))&(process_sorted(i,3)==process_sorted(j,5))
				getindex{i}(count+1)=j;
				count=count+1;
			end
		end
		getindex{i}(1)=count-1;
	end
end