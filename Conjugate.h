/*
   Conjugate gradient method for diagonalizing
   a sparse symmetric matrix
*/
#ifndef CONJUGATE_H
#define CONJUGATE_H
class Conjugate {
public:
	Conjugate(const int &Dim) ;
	~Conjugate() ;
	//void Start(char y_n) ;
	//void SetToZero(double *f) ;

	double ErrorBar ;

	void NormTo1 ( double *f ) ; // normalize <f|f> to 1
 	double f1timesf2(const double *f, const double *g) ;
	int abc_2(const long &iter) ;
	void abc_4() ;

	double test;
	double *f0 ;
	double *f1 ;
	double *f2 ;
	double *f3 ;
	long Dim ;
	double eng ;   // EigenValue of Hamiltonian
	void Ini_f0() ;

private:

	double y00 ;
	double y01 ;
	double y02 ;
	double y22 ;
	double x33_old ;
	inline void CreateSpace() ;
    inline void FreeSpace() ;

} ;
#endif
