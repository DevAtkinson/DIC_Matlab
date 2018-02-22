function script
	for i=1:16
		A(:,:,i)=magic(500);
	end
	save('delete_this_soon.mat','A','-v7.3','-nocompression')
	exampleObject = matfile('delete_this_soon.mat');
	for i=1:100
		tic
		checking(exampleObject);
		time(i)=toc;
	end
	mean(time)

	for i=1:100
		tic
		checking2(exampleObject);
		time(i)=toc;
	end
	mean(time)
end

function checking(obj)
	a=[15 6 15 25 14];
	obj.A(1:5,1:5,:);
end

function checking2(obj)
	a=[15 6 15 25 14];
	obj.A(1:16,1:5,1:5);
end