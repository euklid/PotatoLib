#include "point.h"

Point::Point() {}

Point::Point(unsigned int dimension) :
    m_dim(dimension)
{
    m_coord = std::vector<double>(m_dim,0);
}
double Point::dist(Point const & p1, Point const & p2)
{
    assert(p1.get_dimension() == p2.get_dimension());
    unsigned int dim = p1.get_dimension();
    double res = 0;
    int summand = 0;
    for(unsigned int i = 0; i< dim; i++)
    {
        summand = (p1[i]-p2[i]);
        summand *= summand;
        res += summand;
    }
    return sqrt(res);
}

const double Point::operator[](unsigned int coord) const
{
    assert(coord < m_dim);
    assert(m_coord.size() < coord);
    return m_coord[coord];
}

double Point::operator[](unsigned int coord)
{
    assert(coord < m_dim);
    assert(m_coord.size() < coord);
    return m_coord[coord];
}

unsigned int Point::get_dimension() const
{
    return m_dim;
}
