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
  m_epsL = 1e-4;
  m_tw_omega = 1.2;
}

void SuperGrid::print_parameters() const
{
   printf("SuperGrid parameters left=%i, right=%i, x0=%e, x1=%e, width=%e, transition=%e epsL=%e\n", m_left, m_right, m_x0, m_x1, m_width, m_trans_width,m_epsL);
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
      CHECK_INPUT(false, "The supergrid taper functions at the left and right must be separated. Here x0+width = " << m_x0+m_width << 
		  " and x1-width = " << m_x1-m_width);
    }
    
  }
  else if( m_left )
  {
    if (m_x0+m_width > m_x1 )
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the left must be smaller than the domain. Here x0+width = " << m_x0+m_width << 
		  " and x1 = " << m_x1);
    }
  }    
  else if( m_right )
  {
    if (m_x0 > m_x1-m_width )
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the right must be smaller than the domain. Here x0 = " << m_x0 << 
		  " and x1-width = " << m_x1-m_width );
    }
  }
}

double SuperGrid::dampingCoeff(double x) const
{
  double f=0.;
  if (m_left && x < m_x0+m_width)
// the following makes the damping transition in 0 < const_width <= x <= const_width+trans_width = m_width
// constant damping in 0 <= x <= const_width
    f=sigma( (m_x0 + m_width - x)/m_trans_width); 
  else if (m_right && x > m_x1-m_width)
// the following makes the damping transition in m_x1-m_width < x < m_x1 - const_width < m_x1
// constant damping in m_x1 - const_width <= x <= m_x1
    f=sigma( (x - (m_x1-m_width) )/m_trans_width);
  return f;
}


// used for damping coefficient
double SuperGrid::sigma(double xi) const
{
   double f;
   if (xi<=0.)
      f = 0;
   else if (xi>=1.)
      f = 1.0;
   else
//    f=xi*xi*xi*(10 - 15*xi + 6*xi*xi);
//    f = fmin + (1.-fmin)*xi*xi*xi*(10 - 15*xi + 6*xi*xi);
// C4 function
//    f = fmin + (1.-fmin)* xi*xi*xi*xi*xi*( 
//      126 - 420*xi + 540*xi*xi - 315*xi*xi*xi + 70*xi*xi*xi*xi );
// C5 function
      f =  xi*xi*xi*xi*xi*xi*(
    462-1980*xi+3465*xi*xi-3080*xi*xi*xi+1386*xi*xi*xi*xi-252*xi*xi*xi*xi*xi);
   return f;
}

double SuperGrid::stretching( double x ) const
{
   return 1-(1-m_epsL)*dampingCoeff(x);
}

double SuperGrid::tw_stretching( double x ) const
{
   return 1 + 0.5*sin(m_tw_omega*x);
}
  
void SuperGrid::set_twilight( double omega )
{
   m_tw_omega = omega;
}
