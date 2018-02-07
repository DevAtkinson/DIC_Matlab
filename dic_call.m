function dic_call
	close all
	clc
	the_parallel_pool=gcp;
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	inc=5;
	number_of_parallel_evaluations=16;
	loop_parallel=1;
	current_folder=pwd;
	addpath(strcat(current_folder,'\readimxstuff'));
	
	% save_as='Richard_CTC.mat';
	save_as='FixNR_Sunday2_lucasKanade.mat';
	% save_as='Richard_CTL_41_new_corr_order.mat';
	% image_count=max(size(FileName));


	%define subset size
	subsize=41;
	stepsize=20;

	B=[(1+P2), P3, P1;
		P5, (P6+1), P4;
		0 0 1]
	X=[dx;dy;1];

	% symbolic_warp(B,X)
	% Proc.stepsize=stepsize;
	% Proc.subsize=subsize;
	% Proc.Warp=B;
	% Proc.WarpVec=X;
	if exist(save_as,'file')
		load(save_as);
		Proc.correlated_to=1;
		subpos=Proc.subpos;
		guess_store=Proc.guess;
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=Proc.starting_subset;
		process_order=Proc.process_order;
		valid_subsets=Proc.valid_subsets;
		getIndex=Proc.getIndex;
		FileName=Proc.FileName;
		PathName=Proc.PathName;
		% inc=Proc.inc;
		Proc.inc=inc;
		X=Proc.WarpVec;
		B=Proc.Warp;
		stepsize=Proc.stepsize;
		subsize=Proc.subsize;
		current_image=Proc.correlated_to;
		% Proc.stepsize=stepsize;
		% Proc.subsize=subsize;
		% save_as='Richard_CTC_2_subset_41.mat';
		symbolic_warp(B,X)
	else 
		[FileName,PathName] = uigetfile('*.im7','Select the images','MultiSelect','on');
		Proc.FileName=FileName;
		Proc.PathName=PathName;
		Proc.inc=inc;
		Proc.stepsize=stepsize;
		Proc.subsize=subsize;
		Proc.Warp=B;
		Proc.WarpVec=X;
		symbolic_warp(B,X)
		for i=1:2
			% image_folder=strcat(PathName,FileName(i))
			image_folder = fullfile( PathName , FileName{i} );
			I{i}=readimx(image_folder);
		end
		F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
		mask=makeMask(F_in);
		% figure
		% imagesc(F_in)
		% polygon=impoly();
		% Proc.polygon=polygon;
		% mask=createMask(polygon);
		[subpos,mask_subsets,valid_subsets]=mask2subsets(mask,subsize,stepsize);
		Proc.subpos=subpos;
		Proc.mask_subsets=mask_subsets;
		Proc.valid_subsets=valid_subsets;
		[xguess,yguess,subx,suby]=seedPoints(PathName, FileName, subpos, mask_subsets,stepsize,subsize);
		Proc.guess=[xguess,yguess];
		guess_store=Proc.guess;
		% guess=[guess_store(1),0,0,guess_store(2),0,0];
		guess=[0,0,0,0,0,0];
		starting_subset=whichSubpos(subpos,stepsize,subx,suby);
		Proc.starting_subset=starting_subset;
		[process_order,getIndex]=correlationOrderUpdated(subpos,starting_subset,valid_subsets);
		Proc.process_order=process_order;
		Proc.getIndex=getIndex;
		Proc.correlated_to=1;
		current_image=1;
		save(save_as,'Proc')
	end

	image_folder = fullfile( PathName , FileName{1} );
	I{1}=readimx(image_folder);
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
    [r_F,c_F]=size(F_in);
 	image_count=max(size(FileName));
 	save_as='FixNR_Sunday2_lucasKanade_ZMN.mat';
 	% clear Proc.im

	
	
	% counter_image=1;
	
	% current_image=1;

	elements=sum(sum(valid_subsets))
	size(process_order)

	for k=(current_image+inc):inc:image_count
		fprintf('image %d\n',k);
		tic
		Proc.im{k}.D=process_order;
		if (k==(1+inc))
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});
			[PP(1,:),Corrr(1)]=DICtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(1,2),process_order(1,3)},'guess',guess,'correlation',3);
			for i=2:elements
				[PP(i,:),Corrr(i)]=DICtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(i-1,:),'correlation',3);
			end
			for j=1:elements
				Proc.im{k}.D(j,6:11)=PP(j,:);
				Proc.im{k}.D(j,12)=Corrr(j);
			end
			
		else
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});
			% ppm = ParforProgMon('Correlation progress', elements,1,500,200);
			parfor i=1:elements
				[PP(i,:),Corrr(i)]=DICtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',Proc.im{k-inc}.D(i,6:11),'correlation',3);
				% ppm.increment();
			end
			for j=1:elements
				Proc.im{k}.D(j,6:11)=PP(j,:);
				Proc.im{k}.D(j,12)=Corrr(j);
			end
			
		end

		toc
		Proc.correlated_to=k;
		save(save_as,'Proc');
	end



end