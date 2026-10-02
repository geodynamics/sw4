#ifndef SW4_GMGHDF5HELPERS_H
#define SW4_GMGHDF5HELPERS_H

#include "Require.h"

static inline double gmg_clamp_coordinate(double value, double lower,
                                          double upper, double tolerance,
                                          long long& outside_count) {
  if (value < lower) {
    if (lower - value > tolerance) outside_count++;
    return lower;
  }
  if (value > upper) {
    if (value - upper > tolerance) outside_count++;
    return upper;
  }
  return value;
}

#ifdef USE_HDF5
#include "hdf5.h"

static inline bool read_hdf5_attr_optional_f64(hid_t loc, const char* name,
                                               double& data)
{
  if (H5Aexists(loc, name) <= 0)
    return false;

  hid_t attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  int ierr = H5Aread(attr_id, H5T_IEEE_F64LE, &data);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);
  return true;
}

static inline void read_gmg_surface_spacing(hid_t dataset_id, double& hx,
                                            double& hy)
{
  double h = 0.0;
  const bool has_h =
      read_hdf5_attr_optional_f64(dataset_id, "resolution_horiz", h);
  const bool has_hx =
      read_hdf5_attr_optional_f64(dataset_id, "x_resolution", hx);
  const bool has_hy =
      read_hdf5_attr_optional_f64(dataset_id, "y_resolution", hy);

  CHECK_INPUT(has_h || (has_hx && has_hy),
              "ERROR: GMG surface dataset must define resolution_horiz or "
              "both x_resolution and y_resolution");

  if (!has_hx)
    hx = h;
  if (!has_hy)
    hy = h;

  CHECK_INPUT(hx > 0 && hy > 0,
              "ERROR: GMG surface spacing must be positive, got hx="
                  << hx << " hy=" << hy);
}

enum GMGSurfacePurpose {
  GMG_SURFACE_TOPOGRAPHY,
  GMG_SURFACE_MODEL_TOP
};

static inline hid_t open_gmg_surface_dataset(hid_t group_id,
                                             const char** surface_name,
                                             GMGSurfacePurpose purpose)
{
  const char* topography_candidates[] = {"topography_bathymetry",
                                         "top_surface"};
  const char* model_top_candidates[] = {"top_surface",
                                        "topography_bathymetry"};
  const char** candidates = purpose == GMG_SURFACE_TOPOGRAPHY
                                ? topography_candidates
                                : model_top_candidates;
  const int ncandidates = 2;

  for (int i = 0; i < ncandidates; i++) {
    if (H5Lexists(group_id, candidates[i], H5P_DEFAULT) > 0) {
      if (surface_name)
        *surface_name = candidates[i];
      return H5Dopen(group_id, candidates[i], H5P_DEFAULT);
    }
  }

  return -1;
}

static inline int read_gmg_surface_dims(hid_t dataset_id, hsize_t dims[3])
{
  dims[0] = 0;
  dims[1] = 0;
  dims[2] = 1;

  hid_t dataspace_id = H5Dget_space(dataset_id);
  ASSERT(dataspace_id >= 0);

  const int top_rank = H5Sget_simple_extent_ndims(dataspace_id);
  CHECK_INPUT(top_rank == 2 || top_rank == 3,
              "ERROR: GMG surface dataset must be rank-2 or rank-3, got "
                  << top_rank);

  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
  if (top_rank == 3)
    CHECK_INPUT(dims[2] == 1,
                "ERROR: GMG surface dataset third dimension must be 1, got "
                    << dims[2]);

  H5Sclose(dataspace_id);
  return top_rank;
}

#endif

#endif
