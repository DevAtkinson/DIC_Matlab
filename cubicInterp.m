function val=cubicInterp(A,x,y)
	val=[1, x, x^2, x^3]*A*[1; y; y^2; y^3];

end