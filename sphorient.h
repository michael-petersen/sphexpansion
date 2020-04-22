/*
sphorient.h

functions to read in the orientation files

clean version, MSP 22 April 2020

 */

#include <iostream>
#include <fstream>

using namespace std;


#include "spline.h"

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


void read_orient (string orient_file, SphOrient& orient) {

  ifstream infile;
  infile.open(orient_file, ios::in);
  if (!infile) {
        cout << "Unable to open file";
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

  // construct spline interpolations
  orient.xspline.set_points(orient.time,orient.xcen);
  orient.yspline.set_points(orient.time,orient.ycen);
  orient.zspline.set_points(orient.time,orient.zcen);
  
}

