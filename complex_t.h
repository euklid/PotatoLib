#ifndef COMPLEX_T_H
#define COMPLEX_T_H

#include "assert.h"
#include "point.h"
#include <cmath>

struct complex_t
{
    double real;
    double img;

    complex_t(double real=0, double img=0)
        : real(real), img(img) {}
    
    complex_t(Point const & p)
    {
        assert(p.get_dimension() == 2);
        real = p[0];
        img = p[1];
    }

    complex_t& operator=(const complex_t & a)
    {
        real=a.real;
        img=a.real;
        return *this;
    }

    complex_t operator+(const complex_t & a) const
    {
        return complex_t(a.real+real, a.img+img);
    }

    complex_t& operator+=(const complex_t & a)
    {
        real+=a.real;
        img+=a.real;
        return *this;
    }

    bool operator==(const complex_t& a) const
    {
        return (real == a.real && img == a.img);
    }

    bool operator!=(const complex_t& a) const
    {
        return !(*this == a);
    }

    complex_t operator-(const complex_t & a) const
    {
        return complex_t(real-a.real, img - a.img);
    }

    complex_t& operator-=(const complex_t & a)
    {
        real-=a.real;
        img-=a.img;
        return *this;
    }

    complex_t & operator*=(const complex_t & a)
    {
        double old_real = real;
        real = old_real*a.real + img*a.img;
        img = old_real*a.img + img*a.real;
        return *this;
    }

    complex_t operator*(const complex_t & a) const
    {
        return complex_t(real*a.real-img*a.img, a.img*real + a.real*img);
    }

    complex_t & operator*=(const double & a)
    {
        real*=a;
        img*=a;
        return *this;
    }

    complex_t operator*(const double a) const
    {
        return complex_t(a*real,a*img);
    }

    double abs2() const
    {
        return real*real + img*img;
    }

    double abs() const
    {
        return std::sqrt(real*real+img*img);
    }

    complex_t conj() const
    {
        return complex_t(real,-img);
    }

    complex_t operator/(const double a) const
    {
        assert(a != 0);
        return complex_t(real/a,img/a);
    }

    complex_t & operator/=(const double a)
    {
        assert(a!=0);
        real/=a;
        img/=a;
        return *this;
    }

    complex_t operator/(const complex_t & a) const
    {
        assert(a != complex_t(0,0));
        return ((*this)*a.conj())/a.abs2();
    }

    complex_t & operator/=(const complex_t & a)
    {
        assert(a!= complex_t(0,0));
        double old_real = real;
        double a_abs2 = a.abs2();
        real = old_real*a.real + img*a.img;
        real /= a_abs2;
        img = img * a.real - real*a.img;
        img /= a_abs2;
        return *this;
    }
    
    complex_t & operator-()
    {
        img*=-1;
        real*=-1;
        return *this;
    }
    
    double arg() const
    {
        assert(img != 0 && real !=0);
        return atan2(img, real);
    }
    
    static complex_t log(complex_t const & fst)
    {
        return complex_t(std::log(fst.abs()),fst.arg());
    }
};


#endif
