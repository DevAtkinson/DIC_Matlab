function dic2
	close all
	clc
	the_parallel_pool=gcp;
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	inc=10;
	number_of_parallel_evaluations=16;
	current_folder=pwd;
	addpath(strcat(current_folder,'\readimxstuff'));
	
	save_as='testing5.mat';
	% image_count=max(size(FileName));


	%define subset size
	subsize=61;
	stepsize=29;

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
		subpos=Proc.subpos;
		guess_store=Proc.guess;
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=Proc.starting_subset;
		process_order=Proc.process_order;
		valid_subsets=Proc.valid_subsets;
		getIndex=Proc.getIndex;
		FileName=Proc.FileName;
		PathName=Proc.PathName;
		inc=Proc.inc;
		X=Proc.WarpVec;
		B=Proc.Warp;
		stepsize=Proc.stepsize;
		subsize=Proc.subsize;
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
		figure
		imagesc(F_in)
		polygon=impoly();
		Proc.polygon=polygon;
		mask=createMask(polygon);
		[subpos,mask_subsets,valid_subsets]=mask2subsets(mask,subsize,stepsize);
		Proc.subpos=subpos;
		Proc.mask_subsets=mask_subsets;
		Proc.valid_subsets=valid_subsets;
		save(save_as,'Proc')
		[xguess,yguess,subx,suby]=seedPoints(PathName, FileName, subpos, mask_subsets,stepsize,subsize);
		Proc.guess=[xguess,yguess];
		guess_store=Proc.guess;
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		save(save_as,'Proc')
		starting_subset=whichSubpos(subpos,stepsize,subx,suby);
		Proc.starting_subset=starting_subset;
		save(save_as,'Proc')
		[process_order,getIndex]=correlationOrder(subpos,starting_subset)
		Proc.process_order=process_order;
		Proc.getIndex=getIndex;
		save(save_as,'Proc')
		Proc.correlated_to=1;
	end

	image_folder = fullfile( PathName , FileName{1} );
	I{1}=readimx(image_folder);
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
    [r_F,c_F]=size(F_in);
 	image_count=max(size(FileName));

	
	
	% counter_image=1;
	current_image=21;
	% current_image=Proc.correlated_to;

	elements=sum(sum(valid_subsets))
	size(process_order)
	[gradx,grady]=imgradientxy(F_in,'prewitt');
	for k=(current_image+inc):inc:image_count
		fprintf('image %d',k);
		tic
		Proc.im{k-1}.D=process_order;
		if (k==(1+inc))
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});
			parcall(1)=parfeval(the_parallel_pool,@DICtrackingNew,2,F_in,G_in,subsize,subpos{starting_subset(1),starting_subset(2)},gradx,grady,guess);
			parcount=0;
			parcount2=2;
			parIndex(1)=1;
			while parcount<elements
				% if  ~all([parcall.Read])
					[ind,Params,Cor]=fetchNext(parcall);
					Proc.im{k-1}.D(parIndex(ind),6:11)=Params;
					Proc.im{k-1}.D(parIndex(ind),12)=Cor;
					
					% getIndex{parIndex(ind)};
					% this2=getIndex{parIndex(ind)}(1)
					for i=1:getIndex{parIndex(ind)}(1)
						% this=parIndex(ind)
						% getIndex{parIndex(ind)}
						parIndex(parcount2)=getIndex{parIndex(ind)}(i+1);
						% getIndex{parIndex(ind)}(i+1)
						% D(getIndex{ind}(i+1),2)
						% D(getIndex{ind}(i+1),3)
						parcall(parcount2)=parfeval(the_parallel_pool,@DICtrackingNew,2,F_in,G_in,subsize,subpos{starting_subset(1),starting_subset(2)},gradx,grady,Proc.im{k-1}.D(parIndex(ind),6:11));
						
						parcount2=parcount2+1;
					end
					parcount=parcount+1;
				% end
				% if parcount==920
				% 	Proc.parcall=parcall;
				% 	Proc.parIndex=parIndex;
				% 	save(save_as,'Proc');
				% end
			end
		else
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});

			



			for i=1:number_of_parallel_evaluations
				parcall(i)=parfeval(the_parallel_pool,@DICtrackingNew,2,F_in,G_in,subsize,subpos{starting_subset(1),starting_subset(2)},gradx,grady,Proc.im{k-1-inc}.D(i,6:11));
			% [PP,CCorr]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{out(1),out(2)},'guess',Proc.im{k-1-inc}.D(Index(out(1),out(2)),6:11))
			end
			parcount=0;
			parcount2=number_of_parallel_evaluations+1;
			parIndex(1)=1;
			while parcount<elements-number_of_parallel_evaluations
				if  ~all([parcall.Read])
					[ind,Params,Cor]=fetchNext(parcall);
					Proc.im{k-1}.D(ind,6:11)=Params;
					Proc.im{k-1}.D(ind,12)=Cor;
					% parcount2
					% process_order(parcount2,2)
					% process_order(parcount2,3)
					parcall(parcount2)=parfeval(the_parallel_pool,@DICtrackingNew,2,F_in,G_in,subsize,subpos{starting_subset(1),starting_subset(2)},gradx,grady,Proc.im{k-1-inc}.D(parcount2,6:11));
					parcount2=parcount2+1;
					parcount=parcount+1;
				end
			end
			for i=1:number_of_parallel_evaluations
				[ind,Params,Cor]=fetchNext(parcall);
				Proc.im{k-1}.D(ind,6:11)=Params;
				Proc.im{k-1}.D(ind,12)=Cor;
			end
		end

		toc
		Proc.correlated_to=k;
		save(save_as,'Proc');
	end



end