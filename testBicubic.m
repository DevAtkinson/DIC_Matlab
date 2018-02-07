function testBicubic
	close all
	clc
	% the_parallel_pool=gcp;
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	inc=5;
	number_of_parallel_evaluations=16;
	loop_parallel=1;
	current_folder=pwd;
	addpath(strcat(current_folder,'\readimxstuff'));
	
	% save_as='Richard_CTC.mat';
	save_as='Richard_CTL_41.mat';
	% image_count=max(size(FileName

	[FileName,PathName] = uigetfile('*.im7','Select the images','MultiSelect','on');
	image_count=max(size(FileName));
	% for i=1:2
	% 	% image_folder=strcat(PathName,FileName(i))
	% 	image_folder = fullfile( PathName , FileName{i} );
	% 	I{i}=readimx(image_folder);
	% end
	% F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
	
	Proc.PathName=PathName;
	Proc.FileName=FileName;
	Proc.done=zeros([image_count,1]);
	for i=1:image_count
		fprintf('image %d\n',i);
		if Proc.done(i)==0
			image_folder = fullfile( PathName , FileName{i} );
			G=readimx(image_folder);
			tic
			coef{i}=getBicubicValues(G.Frames{1,1}.Components{1,1}.Planes{1,1});
			Proc.time(i)=toc;
			Proc.coef{i}=coef{i};
			Proc.done(i)=1;
			save('Richard_CTC_coefficients.mat','Proc');
		end
	end
	
	save('bicubicVals.mat','coef')


end