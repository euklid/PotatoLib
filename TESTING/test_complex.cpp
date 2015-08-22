#include "../complex_t.h"
#include <stdlib.h>
#include <iostream>

int main()
{
    complex_t a(3,4), b(2,3), c(5,7);
    
    // summation
    assert(a+b == c);
    complex_t apb = a;
    apb+=b;
    assert(apb == c);
    
    
    complex_t a_(3,-4);
    //conjugate
    assert(a.conj() == a_);
    
    complex_t ab(-6,17);
    //multiplication
    assert(a*b == ab);
    a*=b;
    assert(a == ab);
    
    //absolute2
    a = complex_t(3,4);
    assert(a.abs2() == 25);
    
    //absolute
    assert(a.abs() == 5);
    
    //subtraction
    assert(a-b == complex_t(1,1));
    a-=b;
    assert(a == complex_t(1,1));
    
    //real multiplication
    assert(b*2 == complex_t(4,6));
    b*=2;
    assert(b == complex_t(4,6));
    
    // unary minus
    b = complex_t(1,1);
    c = b;
    b = -b;
    assert(b + c == complex_t(0,0));
    complex_t d = -b +b;
    //std::cout<< d.real << " " << d.img << std::endl;
    assert(d == complex_t(0,0));
    
    // real division
    complex_t e(4,6);
    assert(e/2 == complex_t(2,3));
    e/=2;
    assert(e == complex_t(2,3));
    
    //division
    
    assert(b == (b*a)/a);
    c = b;
    
    std::cout << c.real << " " << c.img << std::endl;
    std::cout << b.real << " " << b.img << std::endl;
    std::cout << a.real << " " << a.img << std::endl;
    b*=a;
    std::cout << b.real << " " << b.img << std::endl;
    b/=a;
    std::cout << b.real << " " << b.img << std::endl;
    assert(b == c);
    
    // arg and log
    std::cout << complex_t(1,1).arg() << std::endl;
    std::cout << complex_t(1,0).arg() << std::endl;
    std::cout << complex_t(-1,1).arg() << std::endl;
    a = complex_t::log(complex_t(-1,1));
    std::cout << a.real <<" " << a.img << std::endl;
    
    
}