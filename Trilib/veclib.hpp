#pragma once

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <values.h>
#include <vector>
#include <algorithm>
#include <array>

typedef std::array<int,2>    Point2I;
typedef std::array<int,3>    Point3I;
typedef std::array<float,2>  Point2F;
typedef std::array<float,3>  Point3F;
typedef std::array<float,4>  Point4F;

typedef std::array<double,2> Point2D;
typedef std::array<double,3> Point3D;
typedef std::array<double,4> Point4D;

typedef std::array<float,3> Array3F;
typedef std::array<float,4> Array4F;

typedef std::array<double,2> Array2D;
typedef std::array<double,3> Array3D;
typedef std::array<double,4> Array4D;

typedef std::array<int,2> Array2I;
typedef std::array<int,3> Array3I;
typedef std::array<int,4> Array4I;

typedef std::array<float,2>  Vec2F;
typedef std::array<float,3>  Vec3F;

typedef std::array<double,2> Vec2D;
typedef std::array<double,3> Vec3D;
typedef std::array<double,4> Vec4D;

namespace JMath
{
template<class T>
inline T max_value( const T &a, const T &b, const T &c)
{
    return std::max(a,std::max(b, c));
}

template<class T>
inline T min_value( const T &a, const T &b, const T &c)
{
    return std::min(a,std::min(b, c));
}

template<class T>
inline T min_value( const T &a, const T &b, const T &c, const T &d)
{
    return std::min(d, std::min(a,std::min(b, c)));
}

template<class T>
inline T length( const std::array<T,3> &A, const std::array<T,3> &B)
{
    double dx = A[0] - B[0];
    double dy = A[1] - B[1];
    double dz = A[2] - B[2];
    return sqrt( dx*dx + dy*dy + dz*dz );
}

template<class T>
inline T length2( const std::array<T,3> &A, const std::array<T,3> &B)
{
    double dx = A[0] - B[0];
    double dy = A[1] - B[1];
    double dz = A[2] - B[2];
    return dx*dx + dy*dy + dz*dz;
}

template<class T>
inline T magnitude( const std::array<T,3> &A )
{
    return sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );
}

template<class T>
inline T dot_product( const std::array<T,3> &A, const std::array<T,3> &B)
{
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template<class T>
inline T dot_product( const std::array<T,2> &A, const std::array<T,2> &B)
{
    return A[0]*B[0] + A[1]*B[1];
}

template<class T>
inline std::array<T,3> cross_product( const std::array<T,3> &A, const std::array<T,3> &B)
{
    std::array<T,3> C;
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
    return C;
}

template<class T>
T random_value(T minVal, T maxVal)
{
    return minVal + drand48()*(maxVal - minVal);
}


template<class T>
inline std::array<T,3> make_vector( const std::array<T,3> &head, const std::array<T,3> &tail)
{
    std::array<T,3> xyz;
    xyz[0] = head[0] - tail[0];
    xyz[1] = head[1] - tail[1];
    xyz[2] = head[2] - tail[2];
    return xyz;
}

///////////////////////////////////////////////////////////////////////////////
template<class T>
inline std::array<T,3> unit_vector( const std::array<T,3> &vec)
{
    double dl  = magnitude(vec);
    std::array<T,3>  uvec;
    uvec[0] = vec[0]/dl;
    uvec[1] = vec[1]/dl;
    uvec[2] = vec[2]/dl;
    return uvec;
}
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline int unit_vector( const Point3D &head, const Point3D &tail)
{
    std::array<T,3> uvec;
    uvec = make_vector( head, tail);
    double dl  = magnitude(uvec);

    uvec[0]  /=  dl;
    uvec[1]  /=  dl;
    uvec[2]  /=  dl;
    return uvec;
}

template<class T>
inline T mean_value( const std::vector<T> &v)
{
    assert( !v.empty() );
    std::vector<T> tmp(v);
    std::sort( tmp.begin(), tmp.end() );
    return tmp[v.size()/2];
}

template<class T>
inline T average_value( const std::vector<T> &v)
{
    size_t nsize = v.size();

    T sum = 0.0;
    for( size_t i = 0; i < nsize; i++)
        sum += v[i];
    T avg = sum/(double)nsize;
    return avg;
}

template<class T>
inline T standard_deviation( const std::vector<T> &v)
{
    assert( !v.empty() );

    size_t nsize = v.size();

    T avg = average_value( v );

    T sum = 0.0;
    for( size_t i = 0; i < nsize; i++) {
        double vm = v[i] - avg;
        sum += vm*vm;
    }
    return sqrt(1.0/(double)(nsize-1)*sum );
}

///////////////////////////////////////////////////////////////////////////////
template<class T>
inline T angle( T x1, T y1, T x2, T y2)
{
    double theta1 = atan2( (double)y1, (double)x1);
    double theta2 = atan2( (double)y2, (double)x2);

    double dtheta = theta2-theta1;

    if( dtheta >  M_PI) dtheta -= 2.0*M_PI;
    if( dtheta < -M_PI) dtheta += 2.0*M_PI;

    return dtheta;
}
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline T angle( const std::array<T,3> &A, const std::array<T,3> &B)
{
    double AB = dot_product(A,B);
    double Am = magnitude(A);
    double Bm = magnitude(B);

    if( Am < 1.0E-15 || Bm < 1.0E-15) return 0.0;

    double x = AB/(Am*Bm);

    if( x > 1.0)  x = 1.0;
    if( x < -1.0) x = -1.0;

    return acos(x);
}

}
