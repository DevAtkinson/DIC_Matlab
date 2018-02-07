function temp
	syms subsize;
	F= sym('F_%d_%d', [subsize subsize]);
	G= sym('G_%d_%d', [subsize subsize]);
	Fmean=mean(mean(F));
	Ftemp=F-ones([subsize, subsize])*Fmean;
	dF_temp=sum(sum(Ftemp.^2));
	dF=sqrt(dF_temp);
	Gmean=mean(mean(G));
	% fprintf('2\n');
	Gtemp=G-ones([subsize, subsize])*Gmean;
	% fprintf('3\n');
	dG_temp=sum(sum(Gtemp.^2));
	% fprintf('4\n');
	dG=sqrt(dG_temp);
	S=sum(sum((Ftemp./dF-Gtemp./dG).^2));
	ans=functionalDerivative(S,G(subsize,subsize))
end