module type_defs
  use iso_c_binding
  integer, parameter:: sp = c_float
  integer, parameter:: dp = c_double
  integer, parameter:: si = c_int
  integer, parameter:: li = c_long
  real(dp), parameter :: pi = acos(-1.d0)
end module type_defs
