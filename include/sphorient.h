/*
sphorient.h

functions to read in the orientation files

MSP 22 Apr 2020 clean version
MSP 25 Apr 2020 add linear interpolation for speed

wishlist:
-handle bad files gracefully

 */

#include <iostream>
#include <fstream>

using namespace std;

// decide on splines
#include "spline.h"
bool splineorient = false;


struct SphOrient
{
  int NUMT;             // the number of timesteps
  vector<double> time;  // the time vector, len NUMT
  vector<double> xcen;  // the xcen vector, len NUMT
  vector<double> ycen;  // the ycen vector, len NUMT
  vector<double> zcen;  // the zcen vector, len NUMT

  tk::spline xspline; // spline representation of xcenter
  tk::spline yspline; // spline representation of ycenter
  tk::spline zspline; // spline representation of zcenter
};



void interpolate_centre(double desired_time,
			SphOrient& orient,
			vector<double>& centre)
{
  // find the nearest time
  // the orient array MUST be evenly spaced, check exp output
  double dt = orient.time[1] - orient.time[0];

  int indx = (int)( (desired_time-orient.time[0])/dt);
  if (indx<0) indx = 0;
  if (indx>orient.NUMT-2) indx = orient.NUMT - 2;

  double x1 = (orient.time[indx+1] - desired_time)/dt;
  double x2 = (desired_time - orient.time[indx])/dt;

  centre[0] = (x1*orient.xcen[indx] + x2*orient.xcen[indx+1]);
  centre[1] = (x1*orient.ycen[indx] + x2*orient.ycen[indx+1]);
  centre[2] = (x1*orient.zcen[indx] + x2*orient.zcen[indx+1]);

}

void spline_centre(double desired_time,
		   SphOrient& orient,
		   vector<double>& centre)
{
  centre[0] = orient.xspline(desired_time);
  centre[1] = orient.yspline(desired_time);
  centre[2] = orient.zspline(desired_time);  
}

void return_centre(double desired_time,
		   SphOrient& orient,
		   vector<double>& centre)
{

  if (splineorient) {
    spline_centre(desired_time, orient, centre);
  } else {
    interpolate_centre(desired_time, orient, centre);
  }
	
}


void read_orient (string orient_file, SphOrient& orient) {

  ifstream infile;
  infile.open(orient_file, ios::in);
  if (!infile) {
        cout << "sphorient::read_orient: Unable to open file";
        exit(1); // terminate with error
  }
  
  string line;
  
  int linenum = 0;

  while (getline(infile, line)) {
    istringstream ss(line);

    if (linenum==0) {
	ss >> orient.NUMT;
	orient.time.resize(orient.NUMT);
	orient.xcen.resize(orient.NUMT);
	orient.ycen.resize(orient.NUMT);
	orient.zcen.resize(orient.NUMT);
      } else {
        ss >> 
        orient.time[linenum-1] >> 
        orient.xcen[linenum-1] >>
        orient.ycen[linenum-1] >>
        orient.zcen[linenum-1];
    }
  linenum ++;

  }

  infile.close();

  // construct spline interpolations if necessary
  if (splineorient) {
  orient.xspline.set_points(orient.time,orient.xcen);
  orient.yspline.set_points(orient.time,orient.ycen);
  orient.zspline.set_points(orient.time,orient.zcen);
  }
  
}

