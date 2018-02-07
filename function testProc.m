function testProc
	a=[1 1 0 1 1;
	1 1 0 1 1;
	1 1 0 0 1;
	1 1 1 1 1;
	1 1 1 1 1];

	[process_sorted,getindex]=correlationOrderUpdate(a,[2,1],a)

end