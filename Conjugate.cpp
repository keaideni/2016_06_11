#include <time.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
using namespace std;
#include "Conjugate.h"

Conjugate::Conjugate(const int &dim) {
	Dim = dim ;
	CreateSpace() ;
	Ini_f0() ;

} ;

/*
void Conjugate::Start(char y_n) {
	CreateSpace() ;
	if ( y_n == 'y' ) Ini_f0() ; // randomly initialize f0 if y_n=='y'
}
*/

void Conjugate::Ini_f0() {
	srand( time(NULL) ) ;
	for ( long i = 0 ; i < Dim ; i++ ){
		f0[i] = (double) rand() ;
	}   
}

int Conjugate::abc_2(const long &iter) {
	y00 = f1timesf2(f0, f0) ;   // y00 = <f0|f0>
	y01 = f1timesf2(f0, f1) ;	// y01 = <f0|H|f0>
	eng = y01 / y00 ;    // energy
	/*
	if(iter%5==0){
		cout << "    Energy[" << iter << "] = " << setprecision(18) << eng << endl ;
	}
	*/
	
	double aux = 2.0 / y00 ;
	for ( long j = 0 ; j < Dim ; j++ )
		f3[j] = ( f1[j] - eng * f0[j] ) * aux ;
	double x33 = f1timesf2(f3, f3) ; // x33=<f3|f3>
	//std::cout<<(y00 * x33) / (eng * eng)<<std::endl;
	if( (y00 * x33) / (eng * eng) < ErrorBar ) return 1 ;
	if( iter == 0 ) {
		for ( long i = 0 ; i < Dim ; i++ )
			f2[i] = - f3[i] ;
	} else {
		double aaa = x33 / x33_old ;
		for ( long k = 0 ; k < Dim ; k++ )
			f2[k] = - f3[k] + aaa * f2[k]  ;
	}
	x33_old = x33 ;
	y02 = f1timesf2( f0, f2 ) ; // y02=<f0|f2>
	y22 = f1timesf2( f2, f2 ) ;  // y22=<f2|f2>
	return 0 ;
}

void Conjugate::abc_4() {
	double x03 = f1timesf2(f0, f3) ;
	double x23 = f1timesf2(f2, f3) ;
	double xa = x23 * y02 - x03 * y22 ;
	double xb = x23 * y00 - y01 * y22 ;
	double xc = x03 * y00 - y01 * y02 ;
	double alpha = (-xb + sqrt( xb * xb - 4.0 * xa * xc)) / (2.0 * xa) ;
	for ( long i = 0 ; i < Dim ; i++ ) {
		f0[i] += alpha * f2[i] ;
		f1[i] += alpha * f3[i] ;
		f3[i] = 0.0 ;
	}
}

double Conjugate::f1timesf2(const double *f, const double *g) {
	double x = 0.0 ;
	for ( long i = 0 ; i < Dim ; i++ )
        x += f[i] * g[i] ;
	return x ;
}




void Conjugate::NormTo1( double *f ) {
	double x = sqrt( f1timesf2(f, f) );
	for ( long i = 0 ; i < Dim ; i++ )
        f[i] /= x ;
}

inline void Conjugate::CreateSpace() {
	f0 = new double[Dim] ; 	assert( f0 != NULL ) ;
	f1 = new double[Dim] ; 	assert( f1 != NULL) ;
	f2 = new double[Dim] ; 	assert( f2 != NULL) ;
	f3 = new double[Dim] ; 	assert( f3 != NULL) ;
}

inline void Conjugate::FreeSpace() {
	delete [] f0 ;
	delete [] f1 ;
	delete [] f2 ;
	delete [] f3 ;
}

Conjugate::~Conjugate() {
// delete [] f0 ;
	FreeSpace() ;
}

