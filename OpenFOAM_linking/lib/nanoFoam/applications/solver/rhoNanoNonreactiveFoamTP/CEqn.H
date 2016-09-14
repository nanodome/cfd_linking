
forAll(mesh.cells(), celli)
{
	//k_tot[celli] = 8.26e04 * Foam::exp(-10681.0 / T[celli]);	      // Titanium (Spicer 2002)
	k_tot[celli] = 1.93E+9 * Foam::exp(-8725.11/T[celli]);		      // Fe(CO)5 => FE + 5CO

	//k_surf[celli] = 4.9e01 * Foam::exp(-8993.0 / T[celli]); 	      //in [m/s]
}

surfaceScalarField phi_f("phi_f",phi / linearInterpolate(rho));


 fvScalarMatrix CEqn
 (
     fvm::ddt(C)
     + fvm::div(phi_f, C)
     - fvm::laplacian(turbulence->muEff() / rho, C)
     + fvm::Sp(k_tot,C)
 );
fvOptions.constrain(CEqn);

CEqn.solve(mesh.solver("C"));
fvOptions.correct(C);

dimensionedScalar N_av("N_av",dimless/dimMoles,6.022140857e23);

C = max(C, 0*C);  	
Info<<"min C"<<min(C)<<endl;
Info<<"max C"<<max(C)<<endl;

concSource = C * k_tot;


