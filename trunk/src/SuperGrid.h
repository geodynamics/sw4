// -*-c++-*-
#ifndef SUPERGRID_H
#define SUPERGRID_H

class SuperGrid 
{

public:
SuperGrid();
void define_taper(bool left, double leftStart, bool right, double rightEnd, 
		  double width );
double sigmaScale(double x) const;
double dampingCoeff(double x) const;
double stretching( double x ) const;
double tw_stretching( double x ) const;
double get_tw_omega() const {return m_tw_omega;}
void   set_twilight( double omega );
void   print_parameters() const;

private:
bool m_left, m_right;
double m_x0, m_x1, m_width, m_trans_width, m_const_width;
double m_epsL, m_tw_omega;
double sigma(double xi) const;

};

#endif
