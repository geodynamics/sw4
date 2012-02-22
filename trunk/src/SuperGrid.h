// -*-c++-*-
#ifndef SUPERGRID_H
#define SUPERGRID_H

class SuperGrid 
{

public:
SuperGrid();
void define_taper(bool left, double leftStart, bool right, double rightEnd, 
		  double width, double transWidth);
double velocityCoeff(double x);
double dampingCoeff(double x);
void print_parameters();

private:
bool m_left, m_right;
double m_x0, m_x1, m_width, m_trans_width, m_const_width;
double phi0(double xi);
double psi0(double xi);

};

#endif
