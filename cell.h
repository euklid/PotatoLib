#ifndef CELL_H
#define CELL_H

#include <vector>
#include <map>
#include "element.h"
#include "point.h"
#include "complex_t.h"
#include <string>

class Cell
{
public:
    Cell();
    Cell(unsigned int dimension, double size, Point const & center);
    virtual const std::vector<double> & get_moments() const;
    virtual const std::vector<complex_t> & get_moments_cmp() const;
    virtual void set_moments(std::vector<double> const & moments);
    virtual void set_moments_cmp(std::vector<complex_t> const & moments_cmp);
    virtual const std::vector<double> & get_local_exps() const;
    virtual const std::vector<complex_t>& get_local_exps_cmp() const;
    virtual void set_local_exps(std::vector<double> const & local_exps);
    virtual void set_local_exps_cmp(std::vector<complex_t> const & local_exps_cmp);
    virtual std::vector<Cell*> & divide() = 0;
    //virtual void set_elements(std::vector<Element*> const & elements);
    virtual const std::vector<Element*> get_elements() const;
    virtual const std::vector<Element*> & get_source_elements() const;
    virtual void set_source_elements(std::vector<Element*> const & el);
    virtual const std::vector<Element*> & get_target_elements() const;
    virtual void set_target_elements(std::vector<Element*> const & el);
	virtual unsigned int number_of_elements() const;
    virtual void set_father(Cell* father);
    virtual Cell * const get_father() const;
    virtual std::vector<Cell*> const & get_children() const ;
    virtual unsigned int get_id() const;
    virtual void set_id(unsigned int index);
    virtual unsigned int get_dimension() const;
    virtual Point const & get_center() const;
    virtual void set_center(Point const & center);
    virtual bool is_leaf() const;
    virtual bool contains_point(Point const & pt) const = 0;
    virtual void set_level(unsigned int lvl);
    virtual unsigned int get_level() const;
    virtual const double get_size() const;
    virtual void add_to_interaction_list(Cell* cell);
    virtual void add_to_direct_list(Cell* cell);
    virtual std::vector<Cell*> get_interaction_list() const;
    virtual std::vector<unsigned int> get_interaction_list_ids() const;
    virtual std::vector<Cell*> get_direct_list() const;
    virtual std::vector<unsigned int> get_direct_list_ids() const;
    virtual bool has_level_grid_position() const;
    virtual Point const & get_level_grid_position() const;
    virtual void set_level_grid_position(Point const & grid_pos);
    
    /**
     * for the preconditioning set diagonal start position of leaf block matrix 
     * in preconditioner. This is equal to the sum of all elements of all leaves
     * that came before this one
     * @param leaf_block_start_pos 0-index based diagonal position of first 
     *        diagonal element of block matrix
     */
    virtual void set_leaf_block_start_pos(int leaf_block_start_pos);
    
    /**
     * see set_leaf_block_start_pos(int leaf_block_start_pos)
     * @return 
     */
    virtual int get_leaf_block_start_pos() const;
    
    virtual int get_leaf_number() const;
    virtual void set_leaf_number(int leaf_number);
    
    virtual std::string debug_info() const;

    
    virtual ~Cell();

protected:
    std::vector<Cell*> m_children;
    Cell* m_father;
    std::map<unsigned int, Cell*> m_interaction_list;
    std::map<unsigned int, Cell*> m_direct_list;
    std::vector<Element*> m_src_elements;
    std::vector<Element*> m_target_elements;
    unsigned int m_id;
    unsigned int m_dim;
    Point m_center;
    Point m_grid_pos;
    std::vector<double> m_moments;
    std::vector<complex_t> m_moments_cmp;
    std::vector<double> m_local_exps;
    std::vector<complex_t> m_local_exps_cmp;
    bool m_has_grid_pos;
    double m_size;
    unsigned int m_level;
    int m_leaf_block_start_pos;
    int m_leaf_number;
};

#endif
