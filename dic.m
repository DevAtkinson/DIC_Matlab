function dic()
	clear
	close all
	clc
	syms dx dy x0 y0 %P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 
	current_folder=pwd;
	addpath(strcat(current_folder,'\readimxstuff'));
	

	save_as='check_passfunc.mat';

	inc=10;
	%define subset size
	subsize=41;
	stepsize=20;

	P = sym('P%d', [6,1])

	B=[(1+P(2)), P(3), P(1);
		P(5), (P(6)+1), P(4);
		0 0 1]
	X=[dx;dy;1];

	correlation_method='NR'

	if exist(save_as,'file')
		load(save_as);
		Proc.correlated_to=1;
		subpos=Proc.subpos;
		% guess_store=Proc.guess;
		guess_store=[0 0];
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=Proc.starting_subset;
		process_order=Proc.process_order;
		valid_subsets=Proc.valid_subsets;
		getIndex=Proc.getIndex;
		FileName=Proc.FileName;
		PathName=Proc.PathName;
		% Proc.stepsize=stepsize;
		stepsize=Proc.stepsize;
		% PathName='G:\Work\Richard data\20160909_Richard_Huchzermeyer_CTC_03\CTC_03\';
		inc=Proc.inc;
		% Proc.inc=inc;
		% X=Proc.WarpVec;
		% B=Proc.Warp;
		stepsize=Proc.stepsize;
		subsize=Proc.subsize;
		current_image=Proc.correlated_to;
		% clear Proc.im
		% Proc.stepsize=stepsize;
		% Proc.subsize=subsize;
		% save_as='Richard_CTC_2_subset_41.mat';
		symbolic_warp(B,X)

		if current_image==1
			[process_order,getIndex]=correlationOrderUpdated(subpos,starting_subset,valid_subsets);
			Proc.process_order=process_order;
		end
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

		xguess=0;
		yguess=0;

		Proc.guess=[xguess,yguess];
		guess_store=Proc.guess;
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=whichSubpos(subpos,stepsize,subx,suby);
		Proc.starting_subset=starting_subset;
		[process_order,getIndex]=correlationOrderUpdated(subpos,starting_subset,valid_subsets);
		Proc.process_order=process_order;
		Proc.getIndex=getIndex;
		Proc.correlated_to=1;
		current_image=1;
		save(save_as,'Proc')
	end

	% save_as='delete_me3.mat';
	% save_as='check_NR_symbolic.mat';

	if strcmp(correlation_method,'NR')==1
		NR_symbolic(B,X,save_as,P);
	elseif strcmp(correlation_method,'LK')==1
		symbolic_warp(B,X);
	end
	% fprintf('DONEEEEEEEEEE\n');
	

	image_folder = fullfile( PathName , FileName{1} );
	I{1}=readimx(image_folder);
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
    [r_F,c_F]=size(F_in);
 	image_count=max(size(FileName));
	
	elements=sum(sum(valid_subsets))
	size(process_order)
	Pend=5+max(size(guess));
	process_order2=process_order(:,2:3);
	[hmm] = unique(process_order2, 'rows', 'first');
	size(hmm)

	% starting=parpool

	%determine name of coefficients saving file based on file save name
	[~,name_short,~] = fileparts(save_as);
	coef_save_name=strcat(name_short,'_coefs.mat');

	for k=(current_image+inc):inc:image_count
		fprintf('image %d\n',k);
		tic
		Proc.im{k}.D=process_order;
		image_folder = fullfile( PathName , FileName{k} );
		I{3}=readimx(image_folder);
		G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});

		% check if the method is newton-raphson and if so then interpolate image
		if strcmp(correlation_method,'NR')==1
			if exist(coef_save_name,'file')
				load(coef_save_name);
				if coef_var.im==k
					coef=coef_var.coef;
				else
					clear coef_var;
					coef=getBicubicValues(G_in);
					coefficient_time=toc
					coef_var.coef=coef;
					coef_var.im=k; 
					save(coef_save_name,'coef_var');
				end
			else
				coef=getBicubicValues(G_in);
				coefficient_time=toc
				coef_var.coef=coef;
				coef_var.im=k; 
				save(coef_save_name,'coef_var');
			end
		end
		
		if (k==(1+inc)) %if this is the first image to be correlated
			tic
			if strcmp(correlation_method,'NR')==1
				[PP(1,:),Corrr(1),iter(1)]=NRtracking3_temp('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'stepsize',stepsize,'subset position',subpos{process_order(1,2),process_order(1,3)},'guess',guess,'coef',coef((subpos{process_order(1,2),process_order(1,3)}.coords(1)-stepsize):(subpos{process_order(1,2),process_order(1,3)}.coords(3)+stepsize),(subpos{process_order(1,2),process_order(1,3)}.coords(2)-stepsize):(subpos{process_order(1,2),process_order(1,3)}.coords(4)+stepsize),:),'coef_shift',[0 0 1],'save_name',save_as)
			elseif strcmp(correlation_method,'LK')==1
				[PP(1,:),Corrr(1)]=DICtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(1,2),process_order(1,3)},'guess',guess,'correlation',3);
			end
			toc
			for i=2:elements
				if strcmp(correlation_method,'NR')==1
					best_guess=bestGuess(Corrr,i-1);
					[PP(i,:),Corrr(i),iter(i)]=NRtracking3_temp('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'stepsize',stepsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(best_guess,:),'coef',coef((subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize),:),'coef_shift',[0 0 1],'save_name',save_as)
				% 	% [PP(i,:),Corrr(i)]=NRtracking3('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(i-1,:),'coef',coef)
				elseif strcmp(correlation_method,'LK')==1
					[PP(i,:),Corrr(i)]=DICtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(i-1,:),'correlation',3);
				end
			end
			for j=1:elements
				Proc.im{k}.D(j,6:Pend)=PP(j,:);
				Proc.im{k}.D(j,Pend+1)=Corrr(j);
			end
		else
			clear PP;
			clear Corrr;
			clear iter;
			% parallel_pool=gcp;
			% ppm = ParforProgMon('Correlation progress', elements,1,500,200);
			if strcmp(correlation_method,'NR')==1
				parfor i=1:elements
					% shift the search area according to the previous displacement estimates
					Pshift=Proc.im{k-inc}.D(i,6:Pend);
					coef_shift=[floor(Pshift(4)), floor(Pshift(1))]; %shift of coefficients (first y then x)
					coef_shift(3)=1;
					% coofs=[(subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize+coef_shift(1)):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize+coef_shift(1)),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize+coef_shift(2)):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize+coef_shift(2))];
					% coef_pass=coef((subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize),:);
					[PP(i,:),Corrr(i),iter(i)]=NRtracking3_temp('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'stepsize',stepsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',Proc.im{k-inc}.D(i,6:Pend),'coef',coef((subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize*coef_shift(3)+coef_shift(1)):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize*coef_shift(3)+coef_shift(1)),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize*coef_shift(3)+coef_shift(2)):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize*coef_shift(3)+coef_shift(2)),:),'coef_shift',coef_shift,'save_name',save_as);
					% ppm.increment();
				end
				PP
				Corrr
				iter
				[out,minindex]=testOutlier(Corrr)
				for i=1:elements
					if out(i)==1
						% shift the search area according to the previous displacement estimates
						Pshift=PP(minindex,:);
						coef_shift=[floor(Pshift(4)), floor(Pshift(1))]; %shift of coefficients (first y then x)
						coef_shift(3)=2;
						% coofs=[(subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize+coef_shift(1)):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize+coef_shift(1)),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize+coef_shift(2)):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize+coef_shift(2))];
						% coef_pass=coef((subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize),:);
						[PP(i,:),Corrr(i),iter(i)]=NRtracking3_temp('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'stepsize',stepsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(minindex,:),'coef',coef((subpos{process_order(i,2),process_order(i,3)}.coords(1)-stepsize*coef_shift(3)+coef_shift(1)):(subpos{process_order(i,2),process_order(i,3)}.coords(3)+stepsize*coef_shift(3)+coef_shift(1)),(subpos{process_order(i,2),process_order(i,3)}.coords(2)-stepsize*coef_shift(3)+coef_shift(2)):(subpos{process_order(i,2),process_order(i,3)}.coords(4)+stepsize*coef_shift(3)+coef_shift(2)),:),'coef_shift',coef_shift,'save_name',save_as);
						% ppm.increment();
					end
				end
				PP
				Corrr
				iter

			elseif strcmp(correlation_method,'LK')==1
				parfor i=1:elements
					[PP(i,:),Corrr(i)]=DICtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',Proc.im{k-inc}.D(i,6:11),'correlation',3);
					% ppm.increment();
				end
			end

			for j=1:elements
				Proc.im{k}.D(j,6:Pend)=PP(j,:);
				Proc.im{k}.D(j,Pend+1)=Corrr(j);
			end
		end
		total_time=toc
		Proc.correlated_to=k;
		% clear ppm;
		save(save_as,'Proc');
	end
end

function out=bestGuess(Cor,n)
	if n>2
		for i=n:-1:2
			if (Cor(i)<0.05)
				out=i;
				break;
			end
		end
		if ~exist('out','var')
			out=n;
			fprintf('The guess used is not ideal\n');
		end
	else
		out=n;
	end
end

%function to determine outliers in correlation value
function [out,minindex]=testOutlier(Cor)
	[minval,minindex]=min(Cor);
	for i=1:max(size(Cor))
		if Cor(i)>4*minval
			out(i)=1;
		else
			out(i)=0;
		end
			
	end
end