#ifndef EW_GEOGRAPHICCOORD_H
#define EW_GEOGRAPHICCOORD_H

#include <string>

class GeographicCoord
{
  friend std::ostream& operator<<(std::ostream& output, const GeographicCoord& g);
 public:
  // -----------------------------------------------------------------
  // class GeographicCoord
  //
  // Provides a container to keep the longitude, latitude, and depth
  // together
  // ------------------------------------------------------------------
  GeographicCoord(double lat, double lon, double depth);
  GeographicCoord();

  double getLatitude() const;
  double getLongitude() const;
  double getDepth() const;

  void resetDepth(double depth) { mDepth = depth; };
 private:
  double mLatitude, mLongitude, mDepth;

};
#endif
