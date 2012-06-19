// -*-c++-*-
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

using namespace std;

class Polynomial{

friend std::ostream& operator<<(std::ostream& output, const Polynomial& s);

public:
Polynomial();
Polynomial(double c[3]);
double coeff(unsigned int q);

// for efficiency and simplicity reasons, we make the coefficients public
double m_c[3];

private:   

};

#endif
