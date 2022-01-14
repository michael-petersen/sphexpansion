/*
cylecoefs.h

functions to handle the preparations for the cylinder coefficient files

MSP 5 May 2020 clean version

now takes raw outcoef files -- need helper functions to check for data corruption from restarts, etc.

todo:
  - convert to eigen arrays

 */
#ifndef CYLCOEFS_H
#define CYLCOEFS_H

#if HAVEYAML
#include "yaml-cpp/yaml.h"	// YAML support
#endif


using namespace std;

// create 2- and 3-d array types
typedef boost::multi_array<double, 3> array_type3;
typedef boost::multi_array<double, 2> array_type2;

struct CylCoefs
{
  int MMAX;            // the number of azimuthal harmonics
  int NORDER;            // the number of radial terms
  int NUMT;            // the number of timesteps

  vector<double> t;    // the time, len NUMT

  array_type3 coscoefs;   // the cosine coefficient table, sized NUMT,MMAX+1,NORDER
  array_type3 sincoefs;   // the   sine coefficient table, sized NUMT,MMAX+1,NORDER

};

void read_coef_file (string& coef_file, CylCoefs& coeftable) {

  ifstream in(coef_file.c_str());

  double tnow;
  bool   newformat = false;

  // read a template version first, then reread with NUMT specified

  // Attempt to read magic number
  //
  unsigned int tmagic;
  unsigned ssize;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  //cout << setw(14) << tmagic << endl;

  if (tmagic == 202004387) {

    newformat = true;

    // YAML size
    //
    in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));


    // Make and read char buffer
    //
    auto buf = std::make_unique<char[]>(ssize+1);
    in.read(buf.get(), ssize);
    buf[ssize] = 0;		// Null terminate


#if HAVEYAML
    YAML::Node node = YAML::Load(buf.get());

    // Get parameters
    //
    coeftable.MMAX    = node["mmax"].as<int>();
    coeftable.NORDER  = node["nmax"].as<int>();
    tnow              = node["time"].as<double>();

#else // compile without libyaml

    std::string yamlblob = buf.get();

    // loop through string parsing
    // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
    // but NOT the top vote getter!
    std::string delim = "\n";
    std::string entry;
    auto start = 0U;
    auto end = yamlblob.find(delim);
    std::string token2,token3;
    while (end != std::string::npos)
    {
        // retreive the entry
        entry = yamlblob.substr(start, end - start);

        // split the entry up to the ':'
        token2 = entry.substr(0, entry.find(":"));

        // split the entry after the ':'
        token3 = entry.substr(entry.find(":")+2,entry.length()-entry.find(":")-2);
        //std::cout<<token2<<"--"<<token3<<endl;

        if (token2.compare("mmax") == 0)   coeftable.MMAX    = std::stoi(token3);
        if (token2.compare("nmax") == 0)   coeftable.NORDER  = std::stoi(token3);
        if (token2.compare("time") == 0)   tnow              = std::stod(token3);

        // advance counters
        start = end + delim.length();
        end = yamlblob.find(delim, start);
    }

    // last one is always 'nmax' which we NEED

    // retreive the entry
    entry = yamlblob.substr(start, end - start);

    // split the entry up to the ':'
    token2 = entry.substr(0, entry.find(":"));

    // split the entry after the ':'
    token3 = entry.substr(entry.find(":")+2,entry.length()-entry.find(":")-2);
    //std::cout<<token2<<"--"<<token3<<endl;

    if (token2.compare("mmax") == 0)   coeftable.MMAX    = std::stoi(token3);
    if (token2.compare("nmax") == 0)   coeftable.NORDER  = std::stoi(token3);

    //std::cout << coeftable.MMAX  << " " << coeftable.NORDER  << std::endl;


#endif


  } else {


    // rewind to the beginning of the file
    //
    in.clear();
    in.seekg(0);

    // first thing in is NUMT,LMAX,NORDER
    in.read((char *)&tnow, sizeof(double));
    in.read((char *)&coeftable.MMAX, sizeof(int));
    in.read((char *)&coeftable.NORDER, sizeof(int));

  }

  int end,tmp,nowpos,bufsize;
  in.seekg (0, ios::end);
  end = in.tellg();

  if (newformat) {
    // guess the number of full timesteps, assuming ssize is always
    // the same (it isn't)
    coeftable.NUMT = end/(((coeftable.MMAX+coeftable.MMAX+1)*coeftable.NORDER)*sizeof(double)  + 2*sizeof(unsigned int) + ssize*sizeof(char));
    bufsize = (((coeftable.MMAX+coeftable.MMAX+1)*coeftable.NORDER)*sizeof(double)  + 2*sizeof(unsigned int) + ssize*sizeof(char));
  } else {
    // compute the number of full timesteps: extra terms are for tnow,MMAX,NORDER ahead of every coefficient set
    coeftable.NUMT = end/(((coeftable.MMAX+coeftable.MMAX+1)*coeftable.NORDER)*sizeof(double)  + sizeof(double) + sizeof(int) + sizeof(int));
    bufsize = (((coeftable.MMAX+coeftable.MMAX+1)*coeftable.NORDER)*sizeof(double)  + sizeof(double) + sizeof(int) + sizeof(int));
  }

  //cout << setw(14) << (coeftable.MMAX*coeftable.NORDER) << setw(14) << sizeof(double) << setw(14) << end << endl;
  //cout << coeftable.NUMT << endl;

  cout << "cylcoefs.read_coef_file: reading NUMT, LMAX, NORDER from file . . . ";
#if FORMAT
  if (newformat) {
    cout << "NEW FORMAT" << endl;
  } else {
    cout << "OLD FORMAT" << endl;
  }
#endif

#if ORDERS
  cout << setw(18) << coeftable.NUMT << setw(18) << coeftable.MMAX << setw(18) << coeftable.NORDER << endl;
#endif

  // reset to the beginning
  in.seekg (0, ios::beg);


  // resize the coefs array appropriately
  coeftable.coscoefs.resize(boost::extents[coeftable.NUMT][coeftable.MMAX+1][coeftable.NORDER]);
  coeftable.sincoefs.resize(boost::extents[coeftable.NUMT][coeftable.MMAX+1][coeftable.NORDER]);
  coeftable.t.resize(coeftable.NUMT);

  // now cycle through each time
  for (int tt=0;tt<coeftable.NUMT;tt++) {

    // check that we aren't too near the end of the file...
    nowpos = in.tellg();
    //cout << setw(14) << nowpos << setw(14) << end << setw(14) << bufsize << endl;
    if (nowpos + bufsize > end) continue;


    if (newformat) {

      // this would be nice to break out as a separate definition...
      in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));
      in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

      //cout << ssize << endl;

      // Make and read char buffer
      //
      auto buf = std::make_unique<char[]>(ssize+1);
      in.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate

#if HAVEYAML
      YAML::Node node = YAML::Load(buf.get());

      // Get parameters
      //
      //cout << node["time"].as<double>() << endl;
      coeftable.t[tt]   = node["time"].as<double>();

#else // compile without libyaml

      std::string yamlblob = buf.get();

      // loop through string parsing
      // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
      // but NOT the top vote getter!
      std::string delim = "\n";
      std::string entry;
      auto start = 0U;
      auto end = yamlblob.find(delim);
      std::string token2,token3;
      while (end != std::string::npos)
      {
          // retreive the entry
          entry = yamlblob.substr(start, end - start);

          // split the entry up to the ':'
          token2 = entry.substr(0, entry.find(":"));

          // split the entry after the ':'
          token3 = entry.substr(entry.find(":")+2,entry.length()-entry.find(":")-2);
          //std::cout<<token2<<"--"<<token3<<endl;

          if (token2.compare("time") == 0)   coeftable.t[tt]   = std::stod(token3);

          // advance counters
          start = end + delim.length();
          end = yamlblob.find(delim, start);
      }

      // last one is always 'nmax' which we can ignore here
      //std::cout<<coeftable.t[tt]<<endl;

#endif

    } else {

      in.read((char *)&coeftable.t[tt], sizeof(double));
      in.read((char *)&tmp, sizeof(int));
      in.read((char *)&tmp, sizeof(int));

    }

    // debug outcoef file
    //cout << "tnow=" << coeftable.t[tt] << " norder=" << tmp << endl;

    for (int m=0; m<=coeftable.MMAX; m++) {

      for (int ir=0; ir<coeftable.NORDER; ir++) {
        in.read((char *)&coeftable.coscoefs[tt][m][ir], sizeof(double));
      }

      if (m) {
	for (int ir=0; ir<coeftable.NORDER; ir++) {
          in.read((char *)&coeftable.coscoefs[tt][m][ir], sizeof(double));
        }
      }

    } // MMAX loop
  } // NUMT loop



  cout << "success!!" << endl;

}



#endif
