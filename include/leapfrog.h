/*
very basic leapfrog integrator to test a circular orbit in the MW

 */

typedef boost::multi_array<double, 2> array_type2;



void mw_leapfrog(SphExpansion* MW,
		 array_type2 mwcoefs,
		 vector<double> xinit,
		 vector<double> vinit,
		 int nint,
		 double dt,
		 array_type2& orbit)
{
  double fx,fy,fz;

  // include the forces for now
  orbit.resize(boost::extents[9][nint]);

  // initialise beginning values
  orbit[0][0] = xinit[0];
  orbit[1][0] = xinit[1];
  orbit[2][0] = xinit[2];
  orbit[3][0] = vinit[0];
  orbit[4][0] = vinit[1];
  orbit[5][0] = vinit[2];

  //now step forward one, using leapfrog (drift-kick-drift) integrator?
  //    https://en.wikipedia.org/wiki/Leapfrog_integration
  //
  int step = 1;
  
  return_forces(MW, mwcoefs,0.0,
		   orbit[0][0],orbit[1][0],orbit[2][0],
		fx,fy,fz);

  orbit[6][0] = fx;
  orbit[7][0] = fy;
  orbit[8][0] = fz;

  int j;

  for (step=1; step<nint; step++) {

    // advance positions
    for (j=0; j<3; j++) {
      orbit[j][step] = orbit[j][step-1]   + (orbit[j+3][step-1]*dt  )  + (0.5*orbit[j+6][step-1]  * (dt*dt));
    }

    // calculate new forces
    return_forces(MW, mwcoefs, 0.0,
		   orbit[0][step],orbit[1][step],orbit[2][step],
		  fx,fy,fz);

    orbit[6][step] = fx;
    orbit[7][step] = fy;
    orbit[8][step] = fz;

    // advance velocities
    for (j=3; j<6; j++) {
      orbit[j][step] = orbit[j][step-1] + (0.5*(orbit[j+3][step-1]+orbit[j+3][step])  * dt );
    }
    
  }
  
}
