// -*-c++-*-
#ifndef FILTER_H
#define FILTER_H

#include <vector>
#include "SecondOrderSection.h"

using namespace std;
// support for derived quantities of the time derivative are not yet implemented
enum FilterType{lowPass, bandPass};

class Filter{

friend std::ostream& operator<<(std::ostream& output, const Filter& s);

public:
Filter(FilterType type, unsigned int numberOfPoles, unsigned int numberOfPasses, double f1, double f2);
~Filter();
void evaluate(int N, double *u, double *mf);
void computeSOS(double dt);
double estimatePrecursor();

FilterType get_type(){return m_type;}
int get_order(){return m_poles;}
int get_passes(){return m_passes;}
double get_corner_freq1(){return m_f1;}
double get_corner_freq2(){return m_f2;}

private:   
Filter();
double realPoleBP(double f1, double f2, double dt, SecondOrderSection *&sos_ptr);
double complexConjugatedPolesBP(double f1, double f2, double dt, double alpha, 
			      SecondOrderSection *&sos1_ptr, SecondOrderSection *&sos2_ptr);
double realPoleLP(double fc, double dt, SecondOrderSection *&sos_ptr);
double complexConjugatedPolesLP(double fc, double dt, double alpha, SecondOrderSection *&sos_ptr);
void a2d(double n[3], double d[3], Polynomial &b, Polynomial &a);

FilterType m_type;
double m_dt, m_f1, m_f2;
unsigned int m_passes, m_poles;
vector<SecondOrderSection*> m_SOSp;
int m_numberOfSOS;
int m_real_poles, m_complex_pairs;
bool m_initialized;
double m_pole_min_re;

};

#endif
