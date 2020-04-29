/*
very basic leapfrog integrator to test a circular orbit in the MW

(just for testing purposes!)

MSP 25 Apr 2020 version with changing coefficients

 */

typedef boost::multi_array<double, 2> array_type2;


void print_orbit(array_type2 orbit,
		 string orbitfile)
{
  ofstream outorbit;
  outorbit.open(orbitfile);

  outorbit << "# t [Gyr]; x [kpc]; y [kpc]; z [kpc]; vx [km/s] ; vy [km/s] ; vz [km/s] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;

  for (int i=0; i<orbit.shape()[1]; i++) {

    outorbit << setw(14) << orbit[9][i];
    
    for (int j=0; j<9; j++) {
      outorbit << setw(14) << orbit[j][i];
    }
    outorbit << endl;
  }

  outorbit.close();
  
}



void leapfrog(SphExpansion* S,
	      array_type2 coefs,
	      vector<double> xinit,
	      vector<double> vinit,
	      int nint,
	      double dt,
	      array_type2& orbit)
{
  /*

   */
  double fx,fy,fz;

  // include the forces for now
  orbit.resize(boost::extents[10][nint]);

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
  
  S->return_forces(S, coefs,
		orbit[0][0],orbit[1][0],orbit[2][0],
		fx,fy,fz);

  orbit[6][0] = fx;
  orbit[7][0] = fy;
  orbit[8][0] = fz;

  int j;

  for (step=1; step<nint; step++) {

    // advance timestep
    orbit[9][step] = dt*step;

    // advance positions
    for (j=0; j<3; j++) {
      orbit[j][step] = orbit[j][step-1]   + (orbit[j+3][step-1]*dt  )  + (0.5*orbit[j+6][step-1]  * (dt*dt));
    }

    // calculate new forces
    S->return_forces(S, coefs,
		  orbit[0][step],orbit[1][step],orbit[2][step],
		  //orbit[0][step],orbit[1][step],0.,
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


void evolving_leapfrog(SphExpansion* S,
	               array_type3 ecoefs,
	               vector<double> xinit,
	               vector<double> vinit,
	               int nint,
	               double dt,
	               array_type2& orbit)
{
  /*

   */
  double fx,fy,fz;

  // include the forces for now
  orbit.resize(boost::extents[10][nint]);

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

  array_type2 tcoefs;
  select_coefficient_time(0., S->coeftable, tcoefs);
  
  S->return_forces(S, tcoefs,
		orbit[0][0],orbit[1][0],orbit[2][0],
		fx,fy,fz);

  orbit[6][0] = fx;
  orbit[7][0] = fy;
  orbit[8][0] = fz;

  int j;

  for (step=1; step<nint; step++) {

    // advance timestep
    orbit[9][step] = dt*step;

    // pull coefficients at this time
    select_coefficient_time(orbit[9][step], S->coeftable, tcoefs);
    
    // advance positions
    for (j=0; j<3; j++) {
      orbit[j][step] = orbit[j][step-1]   + (orbit[j+3][step-1]*dt  )  + (0.5*orbit[j+6][step-1]  * (dt*dt));
    }

    // calculate new forces
    S->return_forces(S, tcoefs,
		  orbit[0][step],orbit[1][step],orbit[2][step],
		  //orbit[0][step],orbit[1][step],0.,
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
