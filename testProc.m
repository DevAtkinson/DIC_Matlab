function testProc
	a=[1 1 0 1 1;
	1 1 0 1 1;
	1 1 0 0 1;
	1 1 1 1 1;
	1 1 1 1 1];

	[process_sorted,getindex]=correlationOrderUpdated(a,[2,1],a)
	size(process_sorted)
	sum(sum(a))

end