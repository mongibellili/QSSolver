model mliqss_test
Real x1, x2;

initial equation
x1 = -1;
x2 = -2;

equation
der(x1) = -20.0*x1 -80.0*x2 + 1600.0;
der(x2) = 1.24*x1 -0.01* x2 + 0.2;
	annotation(

	experiment(
		MMO_Description="",
		MMO_Solver=mLIQSS,
		MMO_PartitionMethod=Metis,
		MMO_Output={x1, x2},
		Jacobian=Dense,
		MMO_BDF_PDepth=1,
		MMO_BDF_Max_Step=0,
		StartTime=0.0,
		StopTime=100,
		Tolerance={1e-2},
		AbsTolerance={1e-5}
	));
end mliqss_test;
