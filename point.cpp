#include "point.h"
#include <cmath>
#include <iosfwd>

Point::Point() {}

Point::Point(unsigned int dimension) :
    m_dim(dimension)
{
    m_coord = std::vector<double>(m_dim,0);
}

Point::Point(const Point& p)
{
    m_dim = p.m_dim;
    m_coord = p.m_coord;
}

double Point::dist(Point const & p1, Point const & p2)
{
    assert(p1.get_dimension() == p2.get_dimension());
    unsigned int dim = p1.get_dimension();
    double res = 0;
    double summand = 0;
    for(unsigned int i = 0; i< dim; i++)
    {
        summand = (p1[i]-p2[i]);
        summand *= summand;
        res += summand;
    }
    return sqrt(res);
}

double Point::max_norm_dist(const Point &p1, const Point & p2)
{
    assert(p1.get_dimension() == p2.get_dimension());
    unsigned int dim = p1.get_dimension();
    double max = 0;
    double diff = 0;
    for (unsigned int i = 0; i< dim; i++)
    {
        diff = p1[i]-p2[i];
        diff = (diff>0)?diff:-diff;
        if (diff  > max)
        {
            max = diff;
        }
    }
    return max;
}

const double & Point::operator[](unsigned int coord) const
{
    assert(coord < m_dim);
    assert(m_coord.size() > coord);
    return m_coord[coord];
}

double & Point::operator[](unsigned int coord)
{
    assert(coord < m_dim);
    assert(m_coord.size() > coord);
    return m_coord[coord];
}

Point Point::operator+(const Point& p1) const
{
    assert(m_dim == p1.m_dim);
    Point res(m_dim);
    for(int i = 0; i<m_dim; i++)
    {
        res.m_coord[i] = m_coord[i] + p1.m_coord[i];
    }
    return res;
}

Point Point::operator-(const Point& p1) const
{
    assert(m_dim == p1.m_dim);
    Point res(m_dim);
    for(int i = 0; i<m_dim; i++)
    {
        res.m_coord[i] = m_coord[i] - p1.m_coord[i];
    }
    return res;
}

bool Point::operator==(const Point &p) const
{
    if (p.m_dim != m_dim) {
        return false;
    }
    for (int i = 0; i<m_dim; i++) {
        if (m_coord[i] != p.m_coord[i]) {
            return false;
        }
    }
    return true;
}

unsigned int Point::get_dimension() const
{
    return m_dim;
}

std::ostream &operator<<( std::ostream &output, const Point &p )
{
    if(p.m_dim == 0) return output;
    output << p.m_coord[0];
    for(int i = 1; i<p.m_dim; i++)
    {
        output << " " << p.m_coord[i];
    }
    return output;
}
