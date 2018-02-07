function diclsq
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
	save_as='lsqtest.mat';
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
		% Proc.correlated_to=1;
		subpos=Proc.subpos;
		guess_store=Proc.guess;
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=Proc.starting_subset;
		process_order=Proc.process_order
		valid_subsets=Proc.valid_subsets;
		getIndex=Proc.getIndex;
		FileName=Proc.FileName;
		PathName=Proc.PathName;
		inc=Proc.inc;
		% Proc.inc=inc;
		X=Proc.WarpVec;
		B=Proc.Warp;
		stepsize=Proc.stepsize;
		subsize=Proc.subsize;
		current_image=Proc.correlated_to;
		% Proc.stepsize=stepsize;
		% Proc.subsize=subsize;
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
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=whichSubpos(subpos,stepsize,subx,suby);
		Proc.starting_subset=starting_subset;
		[process_order,getIndex]=correlationOrder(subpos,starting_subset)
		Proc.process_order=process_order;
		Proc.getIndex=getIndex;
		save(save_as,'Proc')
		Proc.correlated_to=1;
		current_image=1;
	end

	image_folder = fullfile( PathName , FileName{1} );
	I{1}=readimx(image_folder);
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
    [r_F,c_F]=size(F_in);
 	image_count=max(size(FileName));
 	options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','off','FunctionTolerance',1e-40,'StepTolerance',1e-40,'MaxFunctionEvaluations',1e6,'MaxIterations',1000,'UseParallel',false);

	
	% counter_image=1;
	
	% current_image=1;

	elements=sum(sum(valid_subsets))
	size(process_order)

	for k=(current_image+inc):inc:image_count
		fprintf('image %d\n',k);
		tic
		Proc.im{k-1}.D=process_order;
		if (k==(1+inc))
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});
			if loop_parallel==1
				[PP(1,:)]=lsqnonlin(@(P) lsqDIC(F_in,G_in,subpos{process_order(1,2),process_order(1,3)},subsize,P),guess,[],[],options);
				% [PP(1,:),Corrr(1)]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(1,2),process_order(1,3)},'guess',guess);
				for i=2:elements
					[PP(i,:)]=lsqnonlin(@(P) lsqDIC(F_in,G_in,subpos{process_order(i,2),process_order(i,3)},subsize,P),PP(i-1,:),[],[],options);
					% [PP(i,:),Corrr(i)]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(i-1,:));
				end
				for j=1:elements
					Proc.im{k-1}.D(j,6:11)=PP(j,:);
					% Proc.im{k-1}.D(j,12)=Corrr(j);
				end
			else
				parcall(1)=parfeval(the_parallel_pool,@DICtracking,2,'undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{starting_subset(1),starting_subset(2)},'guess',guess);
				parcount=0;
				parcount2=2;
				parIndex(1)=1;
				while parcount<elements
					% if  ~all([parcall.Read])
						[ind,Params,Cor]=fetchNext(parcall);
						Proc.im{k-1}.D(parIndex(ind),6:11)=Params;
						Proc.im{k-1}.D(parIndex(ind),12)=Cor;
						for i=1:getIndex{parIndex(ind)}(1)
							parIndex(parcount2)=getIndex{parIndex(ind)}(i+1);
							parcall(parcount2)=parfeval(the_parallel_pool,@DICtracking,2,'undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(getIndex{parIndex(ind)}(i+1),2),process_order(getIndex{parIndex(ind)}(i+1),3)},'guess',Proc.im{k-1}.D(parIndex(ind),6:11));
							parcount2=parcount2+1;
						end
						parcount=parcount+1;
					% end
				end
			end
		else
			image_folder = fullfile( PathName , FileName{k} );
			I{3}=readimx(image_folder);
			G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});

			if loop_parallel==1
				parfor i=1:elements
					[PP(i,:)]=lsqnonlin(@(P) lsqDIC(F_in,G_in,subpos{process_order(i,2),process_order(i,3)},subsize,P),Proc.im{k-1-inc}.D(i,6:11),[],[],options);
					% [PP(i,:),Corrr(i)]=DICtracking('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',Proc.im{k-1-inc}.D(i,6:11));
				end
				for j=1:elements
					Proc.im{k-1}.D(j,6:11)=PP(j,:);
					% Proc.im{k-1}.D(j,12)=Corrr(j);
				end
			else
				for i=1:number_of_parallel_evaluations
					parcall(i)=parfeval(the_parallel_pool,@DICtracking,2,'undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',Proc.im{k-1-inc}.D(i,6:11));
				end
				parcount=0;
				parcount2=number_of_parallel_evaluations+1;
				parIndex(1)=1;
				while parcount<elements-number_of_parallel_evaluations
					if  ~all([parcall.Read])
						[ind,Params,Cor]=fetchNext(parcall);
						Proc.im{k-1}.D(ind,6:11)=Params;
						Proc.im{k-1}.D(ind,12)=Cor;
						parcall(parcount2)=parfeval(the_parallel_pool,@DICtracking,2,'undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(parcount2,2),process_order(parcount2,3)},'guess',Proc.im{k-1-inc}.D(parcount2,6:11));
						parcount2=parcount2+1;
						parcount=parcount+1;
					end
				end
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