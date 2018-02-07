function [xguess,yguess,subx,suby]=seedPoints(PathName, FileName, subpos, mask,stepsize,subsize)
		for i=1:2
			image_folder = fullfile( PathName , FileName{i} );
			I{i}=readimx(image_folder);
		end
		F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
		G_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
		figure(999)
		subplot(2,2,1);
		imagesc(F_in.*mask);
		colormap(gray)
		hold on;
		% for i=1:max(size(xtick))
		% 	plot([xtick(i) xtick(i)],[ytick(1) ytick(end)],'c')
		% end
		% for i=1:max(size(ytick))
		% 	plot([xtick(1) xtick(end)],[ytick(i) ytick(i)],'c')
		% end
		key=0;
		while(key==0)
			[subx,suby]=ginput(1);
			out_temp=whichSubpos(subpos,stepsize,subx,suby);
			subx=subpos{out_temp(1),out_temp(2)}.coords(2)
			suby=subpos{out_temp(1),out_temp(2)}.coords(1)
			subplot(2,2,2);
			imagesc(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby)-floor(subsize/2):floor(suby)+floor(subsize/2),floor(subx)-floor(subsize/2):floor(subx)+floor(subsize/2)));
			colormap(gray)
			key = waitforbuttonpress;
		end
		subplot(2,2,1);
		plot(subx,suby,'rx')
		subplot(2,2,2);
		imagesc(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(suby)-floor(subsize/2):floor(suby)+floor(subsize/2),floor(subx)-floor(subsize/2):floor(subx)+floor(subsize/2)));
		colormap(gray)
		hold on
		[xf,yf]=ginput(1)
		plot(xf,yf,'rx');
		subplot(2,2,3);
		imagesc(G_in);
		colormap(gray)
		hold on
		key=0;
		while(key==0)
			[xg1,yg1]=ginput(1)
			subplot(2,2,4);
			imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
			colormap(gray)
			key = waitforbuttonpress;
		end
		subplot(2,2,3);
		plot(xg1,yg1,'rx')
		subplot(2,2,4);
		imagesc(I{2}.Frames{1,1}.Components{1,1}.Planes{1,1}(floor(yg1)-floor(subsize/2):floor(yg1)+floor(subsize/2),floor(xg1)-floor(subsize/2):floor(xg1)+floor(subsize/2)));
		colormap(gray)
		[xg,yg]=ginput(1)
		hold on
		plot(xg,yg,'cx')
		xguess=((xf+subx)-(xg+floor(xg1)));
		yguess=((yf+suby)-(yg+floor(yg1)));
end