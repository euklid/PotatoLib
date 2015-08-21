#ifndef POINT_H
#define POINT_H
#include <vector>
#include <assert.h>
#include <math.h>

class Point
{
public:

    Point();

    Point(unsigned int dimension);
    Point(Point const & p);
    static double dist(Point const & p1, Point const & p2);
    static double max_norm_dist(Point const & p1, Point const & p2);
    double & operator[](unsigned int coord);
    const double & operator[](unsigned int coord) const;
    unsigned int get_dimension() const;
    Point operator-(Point const & p1) const;
    Point operator+(Point const & p1) const;
    bool operator==(Point const & p) const;

private:
    std::vector<double> m_coord;
    unsigned int m_dim;
};

#endif
