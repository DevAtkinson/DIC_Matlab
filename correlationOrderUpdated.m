function [process_sorted,getindex]=correlationOrderUpdated(subpos,selected,valid_subsets)
	[r,c]=size(subpos);

	number_of_subsets=sum(sum(valid_subsets));
	done_subsets=zeros(size(subpos));
	process(1,:)=[0,selected(1),selected(2),0,0];
	[per_uniq,done_subsets]=nextPerimeter(selected,done_subsets,valid_subsets);
	for i=2:size(per_uniq,1)+1
		process(i,1)=1;
		process(i,2:end)=per_uniq(i-1,:);
	end
	[~,indexx,~]=unique(process(:,2:3), 'rows', 'first');
	process=process(sort(indexx,1),:);
	count=2;
	while size(process,1)<number_of_subsets
		[per_uniq,done_subsets]=nextPerimeter(per_uniq,done_subsets,valid_subsets);
		num_of_elements=size(process,1);
		for i=(num_of_elements+1):(size(per_uniq,1)+num_of_elements)
			process(i,1)=count;
			process(i,2:5)=per_uniq(i-num_of_elements,:);
		end
		% process
		[~,indexx,~]=unique(process(:,2:3), 'rows', 'first');
		process=process(sort(indexx,1),:);
		count=count+1;
	end



	% % count=zeros(max([r,c]));
	% count=1;
	% for i=1:r
	% 	for j=1:c
	% 		dis=max([abs(selected(1)-i), abs(selected(2)-j)]);
	% 		if subpos{i,j}.include==1
	% 			A(i,j)=dis;
	% 		else
	% 			A(i,j)=inf;
	% 		end
	% 		count=count+1;
	% 	end
	% end

	% count=1;
	% for i=1:r
	% 	for j=1:c
	% 		dis=max([abs(selected(1)-i), abs(selected(2)-j)]);
	% 		if subpos{i,j}.include==1
	% 			if (i~=1)&(A(i-1,j)<A(i,j))
	% 				process(count,:)=[dis,i,j,(i-1),j];
	% 			elseif (j~=1)&(A(i,j-1)<A(i,j))
	% 				process(count,:)=[dis,i,j,i,(j-1)];
	% 			elseif (i~=r)&(A(i+1,j)<A(i,j))
	% 				process(count,:)=[dis,i,j,(i+1),j];
	% 			elseif (j~=c)&(A(i,j+1)<A(i,j))
	% 				process(count,:)=[dis,i,j,i,(j+1)];
	% 			elseif (j~=c)&(i~=r)&(A(i+1,j+1)<A(i,j))
	% 				process(count,:)=[dis,i,j,(i+1),(j+1)];
	% 			elseif (j~=c)&(i~=1)&(A(i-1,j+1)<A(i,j))
	% 				process(count,:)=[dis,i,j,(i-1),(j+1)];
	% 			elseif (j~=1)&(i~=r)&(A(i+1,j-1)<A(i,j))
	% 				process(count,:)=[dis,i,j,(i+1),(j-1)];
	% 			elseif (j~=1)&(i~=1)&(A(i-1,j-1)<A(i,j))
	% 				process(count,:)=[dis,i,j,(i-1),(j-1)];
	% 			end
	% 			count=count+1;
	% 		end
	% 	end
	% end



	process_sorted=sortrows(process);
	process_sorted(1,:)=[0,selected(1),selected(2),0,0];

	for i=1:size(process_sorted,1)
		count=1;
		for j=1:size(process_sorted,1)
			if (process_sorted(i,2)==process_sorted(j,4))&(process_sorted(i,3)==process_sorted(j,5))
				getindex{i}(count+1)=j;
				count=count+1;
			end
		end
		getindex{i}(1)=count-1;
	end
end

function [per_uniq,done_subsets]=nextPerimeter(A,done_subsets,valid_subsets)
	[r,c]=size(A);
	[r_s,c_s]=size(valid_subsets);
	count=1;
	for i=1:r
		per(count,:)=[A(i,1)-1, A(i,2),A(i,1),A(i,2)];
		per(count+1,:)=[A(i,1)+1, A(i,2),A(i,1),A(i,2)];
		per(count+2,:)=[A(i,1), A(i,2)-1,A(i,1),A(i,2)];
		per(count+3,:)=[A(i,1), A(i,2)+1,A(i,1),A(i,2)];
		count=count+4;
	end
	for i=1:size(per,1)
		if (per(i,1)<1)||(per(i,1)>r_s)||(per(i,2)<1)||(per(i,2)>c_s)
			per(i,:)=[0 0 0 0];
		elseif (done_subsets(per(i,1),per(i,2))==1)||(valid_subsets(per(i,1),per(i,2))==0)
			done_subsets(per(i,1),per(i,2));
			valid_subsets(per(i,1),per(i,2));
			per(i,:)=[0 0 0 0];
		end
	end
	percheck=per(:,1:2);
	[~,index,~]=unique(percheck,'rows','stable');
	per_uniq=per(index,:);
	% per_uniq=unique(per,'rows','stable');
	check_zeros=0;
	for i=1:size(per_uniq,1)
		if per_uniq(i,:)==[0 0 0 0]
			for j=i:size(per_uniq,1)-1
				per_uniq(j,:)=per_uniq(j+1,:);
			end
			check_zeros=1;
		end
	end
	if check_zeros==1
		per_uniq=per_uniq(1:end-1,:);
	end
	% per_uniq=unique(per_uniq,'stable');
	for i=1:size(per_uniq,1)
		done_subsets(per_uniq(i,1),per_uniq(i,2))=1;
	end
end