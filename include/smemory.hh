#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <exception>

#define FAST		// Uncommenting this will cause addition and subtraction operations to run slightly faster, but matrices will not be checked for conformability
//#define MATRIXOP	// Leave commented if not using matrix operations +, -, *, +=, -=, *=

class nonconformable: public std::exception
{
	virtual const char* what() const throw()
	{
		return "non-conformable matrices";
	}
} nonconformable;


class zerosize: public std::exception
{
	virtual const char* what() const throw()
	{
		return "cannot allocate array_t with 0 size";
	}
} zerosize;

class BoundsViolation: public std::exception
{
	virtual const char* what() const throw()
	{
		return "violated the bounds of the array";
	}
} BoundsViolation;


class DimensionalViolation: public std::exception
{
	virtual const char* what() const throw()
	{
		return "provided incorrect number of dimensions to subscript operator";
	}
} DimensionalViolation;

template <typename real_t> class array_t
{
	private:
		const unsigned xdim, ydim, zdim, tdim;
		const unsigned size;

		inline void CopyFrom(const array_t &other)
		{
			for(size_t n=0; n!=size; n++)
			{
				array[n]=other.array[n];
			}	
		}

	public:
		
		double * array;
		array_t(const unsigned n1, const unsigned n2 = 1, const unsigned n3 = 1, const unsigned n4 = 1)://, bool isfftw_malloc=0):
		       	xdim{n1}, ydim{n2}, zdim{n3}, tdim{n4},size{n1*n2*n3*n4}
		{
			#if(!defined(FAST))
				if(size<=0) {throw zerosize; return;}
			#endif
			array= new real_t[size];
			//if(isfftw_malloc) array=(double*) fftw_malloc(sizeof(double)*size);
			//else array=new double[size];
		}
				

		void fill(real_t filler)
		{
			for(int n=0; n!=size; n++)
			{
				array[n]=filler;
			}
		}

		array_t(const array_t &other) :
			xdim(other.xdim), ydim(other.ydim), zdim(other.zdim), tdim(other.tdim), size(other.size), array(new real_t[other.size])
		{
			CopyFrom(other);
		}

		array_t& operator=(const array_t &other)
		{
			CopyFrom(other);
			return *this;
		}

		real_t & operator() (const unsigned x)
		{
			#if (!defined(FAST))
				if(dimension>1) throw DimensionalViolation;
				if( x>=xdim) throw BoundsViolation;
			#endif
			return array[x];
		}
		real_t & operator() (const unsigned x, const unsigned y)
		{
			#if (!defined(FAST))
				if(dimension>2) throw DimensionalViolation;
				if( x>=xdim || y>=ydim) throw BoundsViolation;
			#endif
			return array[(y+ydim*x)];
		}
		real_t & operator() (const unsigned x,const unsigned y,const unsigned z) 
		{
			#if (!defined(FAST))
				if(dimension>3) throw DimensionalViolation;
				if( x>=xdim || y>=ydim || z>=zdim) throw BoundsViolation;
			#endif
			return array[(z+zdim*(y+ydim*x))];
		}
		real_t & operator() (unsigned x, unsigned y, unsigned z, unsigned t) 
		{
			#if (!defined(FAST))
				if(dimension>4) throw DimensionalViolation;
				if( x>=xdim||y>=ydim||z>=zdim||t>=tdim) throw BoundsViolation;
			#endif
			return array[t+tdim*(z+zdim*(y+ydim*x))];
		}
		real_t operator() (const unsigned x) const
		{
			#if (!defined(FAST))
				if(dimension>1) throw DimensionalViolation;
				if( x>=xdim) throw BoundsViolation;
			#endif
			return array[x];
		}
		real_t operator() (const unsigned x, const unsigned y) const
		{
			#if (!defined(FAST))
				if(dimension>2) throw DimensionalViolation;
				if( x>=xdim || y>=ydim) throw BoundsViolation;
			#endif
			return array[(y+ydim*x)];
		}
		real_t operator() (const unsigned x,const unsigned y,const unsigned z) const
		{
			#if (!defined(FAST))
				if(dimension>3) throw DimensionalViolation;
				if( x>=xdim || y>=ydim || z>=zdim) throw BoundsViolation;
			#endif
			return array[(z+zdim*(y+ydim*x))];
		}
		real_t operator() (unsigned x, unsigned y, unsigned z, unsigned t) const
		{
			#if (!defined(FAST))
				if(dimension>4) throw DimensionalViolation;
				if( x>=xdim||y>=ydim||z>=zdim||t>=tdim) throw BoundsViolation;
			#endif
			return array[t+tdim*(z+zdim*(y+ydim*x))];
		}
		#if(defined(MATRIXOP))
			template <class real_c>
			friend array_t<real_c> operator+(const array_t<real_c> &a, const array_t<real_c> &b);
			template <class real_c>
			friend array_t<real_c> operator-(const array_t<real_c> &a, const array_t<real_c> &b);
			template<typename Scalar, class real_c>
			friend array_t<real_c> operator*(const array_t<real_c> &a, const Scalar scalar);
			template<typename Scalar, class real_c>	
			friend array_t<real_c> operator*(const Scalar scalar, const array_t<real_c> &a);

			array_t<real_t>& operator+=(const array_t &b)
			{
				#if (!defined(FAST))
					if(size!=b.size || xdim!=b.xdim || ydim!=b.ydim || zdim!=b.zdim || tdim!=b.tdim) throw nonconformable;
				#endif
				for(size_t n; n<size; n++) array[n]+=b.array[n];
				return *this;
			}

			array_t<real_t>& operator-=(const array_t &b)
			{
				#if (!defined(FAST))
					if(size!=b.size || xdim!=b.xdim || ydim!=b.ydim || zdim!=b.zdim || tdim!=b.tdim) throw nonconformable;
				#endif
				for(size_t n; n<size; n++) array[n]+=b.array[n];
				return *this;
			}

			template<typename Scalar>
			array_t& operator*=(const Scalar scalar)
			{
				for(size_t n; n<size; n++) array[n]*=scalar;
				return *this;
			}
		#endif

		~array_t()
		{
			#if(!defined(FAST))
				if(array!=nullptr){delete[] array; array=nullptr;}
			#else
				delete[] array;
			#endif
		}

};
#if(defined(MATRIXOP))
	template <class real_c>
	array_t<real_c> operator+(const array_t<real_c> &a, const array_t<real_c> &b)
	{
		#if (!defined(FAST))
			if(a.size!=b.size || a.xdim!=b.xdim || a.ydim!=b.ydim || a.zdim!=b.zdim || a.tdim!=b.tdim) throw nonconformable;
		#endif

		array_t<real_c> c(a.xdim,a.ydim,a.zdim,a.tdim);
		for(size_t n; n!=a.size; n++) c.array[n]=a.array[n]+b.array[n];
		return c;
	}

	template <class real_c>
	array_t<real_c> operator-(const array_t<real_c> &a, const array_t<real_c> &b)
	{
		#if (!defined(FAST))
		if(a.size!=b.size || a.xdim!=b.xdim || a.ydim!=b.ydim || a.zdim!=b.zdim || a.tdim!=b.tdim) throw nonconformable;
		#endif
		
		array_t<real_c> c(a.xdim,a.ydim,a.zdim,a.tdim);
		for(size_t n; n!=a.size; n++) c.array[n]=a.array[n]-b.array[n];
		return c;
	}

	template<typename Scalar, class real_c>	
	array_t<real_c> operator*(const array_t<real_c> &a, const Scalar scalar)
	{
		array_t<real_c> c(a.xdim,a.ydim,a.zdim,a.tdim);
		for(size_t n; n<c.size; n++) c.array[n]=a.array[n]*scalar;
		return c;
	}

	template<typename Scalar, class real_c>
	array_t<real_c> operator*(const Scalar scalar, const array_t<real_c> &b)
	{
		array_t<real_c> c(b.xdim,b.ydim,b.zdim,b.tdim);
		for(size_t n; n<c.size; n++) c.array[n]=b.array[n]*scalar;
		return c;
	}
#endif



void *smalloc(int n, const char *name){

  	if (n == 0) return NULL;
  		void *ptr = malloc(n);
  	if (ptr == NULL) {
    		char str[128];
    		printf("Failed to allocate %d bytes for array %s",n,name);
  	}
  	return ptr;
}



void sfree(void *ptr){

  	if (ptr == NULL) return;
  	free(ptr);
}

double *create_1d_double_array(int n1, const char *name){

  	double *array = (double *) smalloc(n1*sizeof(double),name);
  	return array;
}



void destroy_1d_double_array(double *array){

  	if (array == NULL) return;
  	sfree(array);
}


double **create_2d_double_array(int n1, int n2, const char *name){

  	double *data = (double *) smalloc(n1*n2*sizeof(double),name);
  	double **array = (double **) smalloc(n1*sizeof(double *),name);

  	int n = 0;
  	for (int i = 0; i < n1; i++){
    		array[i] = &data[n];
    		n += n2;
  	}
  
  	return array;
}



void destroy_2d_double_array(double **array){

  	if (array == NULL) return;
  	sfree(array[0]);
  	sfree(array);
}



double ***create_3d_double_array(int n1, int n2, int n3, const char *name){

  	int i,j;
  
  	double *data = (double *) smalloc(n1*n2*n3*sizeof(double),name);
  	double **plane = (double **) smalloc(n1*n2*sizeof(double *),name);
  	double ***array = (double ***) smalloc(n1*sizeof(double **),name);

  	int n = 0;
  	for (i = 0; i < n1; i++) {
    		array[i] = &plane[i*n2];
    		for (j = 0; j < n2; j++) {
      			plane[i*n2+j] = &data[n];
      			n += n3;
    		}
  	}

  	return array;
}



void destroy_3d_double_array(double ***array){

  	if (array == NULL) return;
  	sfree(array[0][0]);
  	sfree(array[0]);
  	sfree(array);
}



double ****create_4d_double_array(int n1, int n2, int n3, int n4, const char *name){

  	int i,j,k;

  	double *data = (double *) smalloc(n1*n2*n3*n4*sizeof(double),name);
  	double **cube = (double **) smalloc(n1*n2*n3*sizeof(double *),name);
  	double ***plane = (double ***) smalloc(n1*n2*sizeof(double **),name);
  	double ****array = (double ****) smalloc(n1*sizeof(double ***),name);

  	int n = 0;
  	for (i = 0; i < n1; i++) {
    		array[i] = &plane[i*n2];
    		for (j = 0; j < n2; j++) {
      			plane[i*n2+j] = &cube[i*n2*n3+j*n3];
      			for (k = 0; k < n3; k++) {
				cube[i*n2*n3+j*n3+k] = &data[n];
				n += n4;
      			}
    		}
  	}
  	return array;
}



void destroy_4d_double_array(double ****array)
{
  	if (array == NULL) return;
  	sfree(array[0][0][0]);
  	sfree(array[0][0]);
  	sfree(array[0]);
  	sfree(array);
}




double *****create_5d_double_array(int n1, int n2, int n3, int n4, int n5, const char *name){

  int i,j,k,l;
  
  double *data = (double *) smalloc(n1*n2*n3*n4*n5*sizeof(double),name);
  double **cube = (double **) smalloc(n1*n2*n3*n4*sizeof(double *),name);
  double ***plane = (double ***) smalloc(n1*n2*n3*sizeof(double **),name);
  double ****boss = (double ****) smalloc(n1*n2*sizeof(double ***),name);
  double *****array = (double *****) smalloc(n1*sizeof(double ****),name);
  
  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &boss[i*n2];
    
    for (j = 0; j < n2; j++) {
      boss[i*n2+j] = &plane[(i*n2+j)*n3];
      
      for (k = 0; k < n3; k++) {
	plane[(i*n2+j)*n3+k] = &cube[(i*n2*n3+j*n3+k)*n4];
	
	for(l = 0; n < n4 ; l++){
	  cube[(i*n2*n3+j*n3+k)*n4+l] = &data[n];
	  n += n5;
	}
	
      }
    }
  }
  return array;
}


void destroy_5d_double_array(double *****array)
{
  	if (array == NULL) return;
	sfree(array[0][0][0][0]);
  	sfree(array[0][0][0][0]);
  	sfree(array[0][0][0]);
  	sfree(array[0][0]);
	sfree(array[0]);
  	sfree(array);
}

