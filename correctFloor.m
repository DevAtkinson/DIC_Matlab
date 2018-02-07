function correctFloor(fileID)
	fid = fopen(fileID);
	tline = fgetl(fid);
	count=0;
	while ischar(tline)
		count=count+1;
	    lines{count}=tline;
	    tline = fgetl(fid);
	end
	old='D(D(D(D(floor))))';
	new='floor';
	for i=1:count
		lines{i} = replace(lines{i},old,new);
	end
	old='D(D(D(floor)))';
	new='floor';
	for i=1:count
		lines{i} = replace(lines{i},old,new);
	end
	old='D(D(floor))';
	new='floor';
	for i=1:count
		lines{i} = replace(lines{i},old,new);
	end
	old='D(floor)';
	new='floor';
	for i=1:count
		lines{i} = replace(lines{i},old,new);
	end
	fid2=fopen(fileID,'w');
	for i=1:count
		fprintf(fid2,'%s\n',lines{i});
	end
end 