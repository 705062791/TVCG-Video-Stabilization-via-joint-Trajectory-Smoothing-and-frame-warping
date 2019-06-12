/*
	=================================
	Tuple          2009/11/11
	=================================

	There are many kinds of arrays such as array_1d,array_2d,array_3d.
	The elements of these arrays can be a tuple.
	A tuple' elements can be any type and the number of its elements
	can be arbitrary.
*/
#pragma  once

#include <cmath>

#define ForIndex(I,N)							for(int I = 0;      I< int(N); I++)
#define ForRange(I,A,B)							for(int I = (int)A; I<=int(B); I++)
#define ForRangeReverse(I,A,B)                  for(int I = (int)A; I>=int(B); I--)
#define ForIndex(I,N)							for(int I = 0;      I< int(N); I++)

namespace LibIV
{
	namespace Math
	{
		template<typename T_Type, int T_N>
		class Tuple
		{
		protected:
			T_Type m_Array[T_N];

		public:
			Tuple()			{}
			~Tuple()        {}
			Tuple(const T_Type& v)
			{
				//iv_dbg_assert((m_Array!=NULL));
				ForIndex(i,T_N)
					m_Array[i] = v;
			}

			Tuple(const Tuple<T_Type,T_N> & tup)
			{
				//iv_dbg_assert((m_Array!=NULL));
				ForIndex(i,T_N)
					m_Array[i] = tup[i];
			}

			// ===============================================
			// Accessors
			const T_Type& operator[] ( int i ) const
			{
				//iv_dbg_assert((i>=0) && (i<T_N));
				return m_Array[i];
			}

			T_Type & operator[] (int i)
			{
				//iv_dbg_assert((i>=0) && (i<T_N));
				return m_Array[i];
			}
			// ==============================================
			// Assignment
			Tuple<T_Type,T_N>& operator= (const Tuple<T_Type,T_N>& tup)
			{
				//iv_dbg_assert((m_Array!=NULL));
				ForIndex(i,T_N)
					m_Array[i] = tup[i];
				return (*this);
			}

			Tuple<T_Type,T_N>& operator= (const T_Type& v)
			{
				//iv_dbg_assert((m_Array!=NULL));
				ForIndex(i,T_N)
					m_Array[i] = v;
				return (*this);
			}

			bool operator ==(const Tuple<T_Type,T_N> & tup)
			{
				ForIndex(i,T_N)
					if(tup[i] != m_Array[i])
						return false;
				return true;
			}

			// ==============================================
			// Add,subtract,multiply,divide
			Tuple<T_Type,T_N> operator+(const Tuple<T_Type,T_N> & rhs)
			{
				Tuple<T_Type,T_N> ans;
				ForIndex(i,T_N)
					ans[i] = m_Array[i] + rhs[i];
				return ans;
			}

			Tuple<T_Type,T_N> operator+(const Tuple<T_Type,T_N> & rhs) const
			{
				Tuple<T_Type,T_N> ans;
				ForIndex(i,T_N)
					ans[i] = m_Array[i] + rhs[i];
				return ans;
			}

			Tuple<T_Type,T_N> operator-(const Tuple<T_Type,T_N> & rhs)
			{
				Tuple<T_Type,T_N> ans;
				ForIndex(i,T_N)
					ans[i] = m_Array[i] - rhs[i];
				return ans;
			}

			Tuple<T_Type,T_N> operator-(const Tuple<T_Type,T_N> & rhs) const
			{
				Tuple<T_Type,T_N> ans;
				ForIndex(i,T_N)
					ans[i] = m_Array[i] - rhs[i];
				return ans;
			}

			Tuple<T_Type,T_N> operator*(const double rhs)
			{
				Tuple<T_Type,T_N> ans;

				ForIndex(i,T_N)
					ans[i] = m_Array[i] * rhs;

				return ans;
			}

			Tuple<T_Type,T_N> operator*(const double rhs) const
			{
				Tuple<T_Type,T_N> ans;
				ForIndex(i,T_N)
					ans[i] = m_Array[i] * rhs;

				return ans;
			}

			const T_Type *	raw()			const	{return m_Array;}

			T_Type *		raw()					{return m_Array;}

		};//End of Tuple

		typedef Tuple<int,2>		v2i;
		typedef Tuple<int,3>		v3i;
		typedef Tuple<int,4>        v4i;
		typedef Tuple<int,5>        v5i;

		typedef Tuple<float,2>		v2f;
		typedef Tuple<float,3>      v3f;


		typedef Tuple<double,2>		v2d;
		typedef Tuple<double,3>		v3d;
		typedef Tuple<double,4>     v4d;
		typedef Tuple<double,6>     v6d;
		typedef Tuple<double,10>    v10d;
		typedef Tuple<double,16>    v16d;

		typedef Tuple<char,2>		v2c;
		typedef Tuple<char,3>		v3c;

		typedef Tuple<unsigned char,2>		v2b;
		typedef Tuple<unsigned char,3>		v3b;


		inline v2i _v2i_(int x,int y)
		{
			v2i v;
			v[0] = x;
			v[1] = y;
			return v;
		}

		inline v3i _v3i_(int x,int y,int z)
		{
			v3i v;
			v[0] = x;
			v[1] = y;
			v[2] = z;
			return v;
		}

		inline v4i _v4i_(int a ,int b, int c, int d)
		{
			v4i v;
			v[0] = a;
			v[1] = b;
			v[2] = c;
			v[3] = d;
			return v;
		}

		inline v2f _v2f_(float x,float y)
		{
			v2f v;
			v[0] = x;
			v[1] = y;
			return v;
		}

		inline v3f _v3f_(float x,float y,float z)
		{
			v3f v;
			v[0] = x;
			v[1] = y;
			v[2] = z;
			return v;
		}

		inline v2d _v2d_(double x,double y)
		{
			v2d v;
			v[0] = x;
			v[1] = y;
			return v;
		}

		inline v3d _v3d_(double a, double b, double c)
		{
			v3d v;
			v[0] = a;
			v[1] = b;
			v[2] = c;
			return v;
		}

		inline v4d _v4d_(double a, double b, double c, double d)
		{
			v4d v;
			v[0] = a;
			v[1] = b;
			v[2] = c;
			v[3] = d;
			return v;
		}

		inline v6d _v6d_(double a, double b, double c, double d, double e, double f)
		{
			v6d v;
			v[0] = a;
			v[1] = b;
			v[2] = c;
			v[3] = d;
			v[4] = e;
			v[5] = f;
			return v;
		}

		inline v10d _v10d_(double a, double b, double c, double d, double e, double f,
                   double g, double h, double i, double j)
		{
			v10d v;
			v[0] = a;
			v[1] = b;
			v[2] = c;
			v[3] = d;
			v[4] = e;
			v[5] = f;
			v[6] = g;
			v[7] = h;
			v[8] = i;
            v[9] = j;
			return v;
		}


		inline bool isEqualV3i(const v3i& a, const v3i& b)
		{
			return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
		}

		inline void crossProduct_v3d(v3d a ,v3d b, v3d & c)
		{
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];
		}

		inline v3d crossProduct_v3d(v3d a, v3d b)
		{
			v3d c;
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];
			return c;
		}

		inline double dotProduct_v3d(v3d a, v3d b)
		{
			return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
		}

		inline v3d normalize_v3d(const v3d & a)
		{
			double lg = sqrt(dotProduct_v3d(a,a));
			v3d ans;
			ans[0] = a[0]/lg;
			ans[1] = a[1]/lg;
			ans[2] = a[2]/lg;
			return ans;
		}

		inline double square_v3d(const v3d & a)
		{
			return dotProduct_v3d(a,a);
		}

		inline double square_v2d(const v2d & a)
		{
			return (a[0] * a[0] + a[1] * a[1]);
		}

		inline double dist_v2d(const v2d & a, const v2d & b)
		{
			return sqrt(square_v2d(_v2d_(a[0]-b[0],a[1]-b[1])));
		}
	} // End of Math
} // End of LibIV

using namespace LibIV::Math;