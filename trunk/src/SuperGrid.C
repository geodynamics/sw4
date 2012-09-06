#include "SuperGrid.h"
#include "Require.h"
#include <cstdio>
using namespace std;

SuperGrid::SuperGrid()
{
  m_left = false;
  m_right = false;
  m_x0=0.;
  m_x1=1.;
  m_width=0.1;
  m_const_width=0.;
}

void SuperGrid::print_parameters()
{
  printf("SuperGrid parameters left=%i, right=%i, x0=%e, x1=%e, width=%e, transition=%e\n", m_left, m_right, m_x0, m_x1, m_width, m_trans_width);
}

void SuperGrid::define_taper(bool left, double leftStart, bool right, double rightEnd, double width, double transWidth)
{
  m_left = left;
  m_x0 = leftStart;
  m_right = right;
  m_x1 = rightEnd;
  m_width = width;
  m_trans_width = transWidth;
  m_const_width = m_width - m_trans_width;
  
// sanity checks
  if (m_left || m_right)
  {
     double dlen = m_x1-m_x0;
     CHECK_INPUT(m_width > 0., "The supergrid taper width must be positive, not = " << m_width);
     CHECK_INPUT(m_width < dlen, "The supergrid taper width must be smaller than the domain, not = " << m_width);
     CHECK_INPUT(m_trans_width > 0., "The supergrid taper transition width must be positive, not = " << m_trans_width);
     CHECK_INPUT(m_const_width >= 0., "The supergrid const_width = width - trans_width must be non-negative, not = " << m_const_width);
  }
  
  if (m_left && m_right)
  {
    if (m_x0+m_width > m_x1-m_width)
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the left and right must be separated. Here x0+widht = " << m_x0+m_width << 
		  " and x1-width = " << m_x1-m_width);
    }
    
  }
  
}

double SuperGrid::velocityCoeff(double x)
{
  double f=1.;
  if (m_left)
// using trans_width here reduces the velocity in 0<= x <= m_trans_width < m_width
//    f*=phi0( (x - (m_x0))/m_trans_width);
    f*=phi0( (x - (m_x0+m_const_width))/m_trans_width); 
  if (m_right)
// using trans_width here reduces the velocity in m_x1 - m_trans_width <= x <= m_x1
//    f*=phi0( ((m_x1) - x)/m_trans_width);
    f*=phi0( ((m_x1-m_const_width) - x)/m_trans_width);
  return f;
}

double SuperGrid::dampingCoeff(double x)
{
  double f=1.;
  if (m_left)
// the following makes the damping transition in 0 < const_width <= x <= const_width+trans_width = m_width
// constant damping in 0 <= x <= const_width
    f*=psi0( (x - (m_x0+m_const_width))/m_trans_width); 
  if (m_right)
// the following makes the damping transition in m_x1-m_width < x < m_x1 - const_width < m_x1
// constant damping in m_x1 - const_width <= x <= m_x1
    f*=psi0( ((m_x1-m_const_width) - x)/m_trans_width);
  return 1.-f;
}

// used for scaling down velocities
double SuperGrid::phi0(double xi)
{
  double f, fmin=1.e-4;
  if (xi<=0.)
//     f = 0.; // the free surface boundary condition does not like the velocity to be zero on the free surface
    f = fmin;
  else if (xi>=1.)
    f = 1.0;
  else
//    f=xi*xi*xi*(10 - 15*xi + 6*xi*xi);
//    f = fmin + (1.-fmin)*xi*xi*xi*(10 - 15*xi + 6*xi*xi);
     f = fmin + (1.-fmin)* xi*xi*xi*xi*xi*( 
       126 - 420*xi + 540*xi*xi - 315*xi*xi*xi + 70*xi*xi*xi*xi );
  

  return f;
}

// used for damping coefficient
double SuperGrid::psi0(double xi)
{
  double f, fmin=0.;
  if (xi<=0.)
    f = fmin;
  else if (xi>=1.)
    f = 1.0;
  else
//    f=xi*xi*xi*(10 - 15*xi + 6*xi*xi);
//    f = fmin + (1.-fmin)*xi*xi*xi*(10 - 15*xi + 6*xi*xi);
    f = fmin + (1.-fmin)* xi*xi*xi*xi*xi*( 
      126 - 420*xi + 540*xi*xi - 315*xi*xi*xi + 70*xi*xi*xi*xi );
  return f;
}


  
