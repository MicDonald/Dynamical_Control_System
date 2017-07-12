class Interpolation {
	const double* data;
	double spacing;
	unsigned nMax;

public:
	Interpolation (
		const double* data,
		double spacing,
		unsigned nMax
	) noexcept :
		data( data ),
		spacing( spacing ),
		nMax( nMax )
	{}


	double
	interpolate ( // Newton interpolation
		double x
	) const
	{
		unsigned i = x/spacing;
		unsigned id1 = i-1;
		unsigned id2 = i;
		unsigned id3 = i+1;
		unsigned id4 = i+2;
		if(nMax <= id4) throw 0.;

		const double& b0 = data[id1];
		double d1 = (data[id2]-data[id1])/spacing;
		double d2 = (data[id3]-data[id2])/spacing;
		double d3 = (data[id4]-data[id3])/spacing;
		double dSq1 = (d2-d1) / (spacing*2.);
		double dSq2 = (d3-d2) / (spacing*2.);
		double dCu1 = (dSq2-dSq1) / (spacing*3.);

		return b0 + d1*(x-spacing*id1) + dSq1*(x-spacing*id1)*(x-spacing*id2) + dCu1*(x-spacing*id1)*(x-spacing*id2)*(x-spacing*id3);
	}


	double
	diff (	double x ) const
	{
		unsigned i = x/spacing;
		unsigned id1 = i-1;
		unsigned id2 = i;
		unsigned id3 = i+1;
		unsigned id4 = i+2;
		if(nMax <= id4) throw 0.;

		double d1 = (data[id2]-data[id1])/spacing;
		double d2 = (data[id3]-data[id2])/spacing;
		double d3 = (data[id4]-data[id3])/spacing;
		double dSq1 = (d2-d1) / (spacing*2.);
		double dSq2 = (d3-d2) / (spacing*2.);
		double dCu1 = (dSq2-dSq1) / (spacing*3.);

		return d1 + dSq1*(2.*x-(id1+id2)*spacing) + dCu1*( (x-spacing*id2)*(x-spacing*id3) + (x-spacing*id1)*(x-spacing*id3) + (x-spacing*id1)*(x-spacing*id2) );
	}


	double
	diffSq (
		double x
	) const
	{
		unsigned i = x/spacing;
		unsigned id1 = i-1;
		unsigned id2 = i;
		unsigned id3 = i+1;
		unsigned id4 = i+2;
		if(nMax <= id4) throw 0.;

		double d1 = (data[id2]-data[id1])/spacing;
		double d2 = (data[id3]-data[id2])/spacing;
		double d3 = (data[id4]-data[id3])/spacing;
		double dSq1 = (d2-d1) / (spacing*2.);
		double dSq2 = (d3-d2) / (spacing*2.);
		double dCu1 = (dSq2-dSq1) / (spacing*3.);

		return dSq1*2. + 2.*dCu1*( (x-spacing*id1) + (x-spacing*id2) + (x-spacing*id3) );
	}
};

