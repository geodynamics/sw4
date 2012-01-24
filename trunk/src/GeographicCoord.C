#include "Require.h"
#include "GeographicCoord.h"

#include <iostream>
using namespace std;

GeographicCoord::
GeographicCoord(double lat,
		   double lon,
		   double depth):
  mLatitude(lat),
  mLongitude(lon),
  mDepth(depth){

  REQUIRE2(lat >= -90.0 && lat <= 90.0,
	   "latitude must be between -90 to 90 degrees, not" << lat);
  REQUIRE2(lon >= -180.0 && lon <= 180.0,
	   "longitude must be between -180 to 180 degrees, not " << lon);
}

GeographicCoord::
GeographicCoord():
  mLatitude(0.0),
  mLongitude(0.0),
  mDepth(0.0)
{
}

double
GeographicCoord::
getLatitude() const 
{
  return mLatitude;
}

double
GeographicCoord::
getLongitude() const 
{
  return mLongitude;
}

double
GeographicCoord::
getDepth() const 
{
  return mDepth;
}

std::ostream& operator<<(std::ostream& output, const GeographicCoord& g) {
  return output << "(" << g.getLatitude() << "'," << g.getLongitude() << "') @ " << g.getDepth() << " m";
}
