#include <iostream>
#include "tree.h"
#include "element.h"
#include "point_element.h"
#include "point.h"
#include "fmm2d.h"
#include "kernel_laplace_point_2d.h"
#include "fmm_gmres_solver.h"


void read_fmm_config(char* filename,
                     unsigned int & exp_terms,
                     unsigned int & loc_terms,
                     unsigned int & max_cell_elements)
{
    std::ifstream file;
    file.open(filename);
    assert(file.is_open());
    file >> exp_terms >> loc_terms >> max_cell_elements;
    file.close();
}

void read_point2d_elements(char* filename, 
                           std::vector<Element*> & src_el, 
                           std::vector<Element*> & tgt_el,
                           bool & src_eq_tgt)
{
    std::ifstream file;
    file.open(filename);
    assert(file.is_open());
    /**
     * configuration file is of form:
     * 
     * #src_eq_tgt 
     * 0
     * # x(double) y(double) s(0|1) t(0|1) init_val
     * 0.2 0.4 1 0 0.03
     * 0.4 0.3 0 1 0
     * 
     */
    
    file >> src_eq_tgt;
    double x, y, init_val;
    int source, target;
    int id = 0;
    while(file >> x >> y >> source >> target >> init_val)
    {
        Point position(2);
        position[0] = x;
        position[1] = y;
        int type = source*Element::SOURCE + target*Element::TARGET;
        PointElement* el = new PointElement(position,id,type);
        if(source) 
        {
            src_el.push_back(el);
            el->set_value(init_val);
        }
        if(target)
        {
            tgt_el.push_back(el);
        }
        id++;
    }
    file.close();
    
}

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        std::cerr << "no element configuration provided" << std::endl;
        return -1;
    }
    
    if(argc < 3) 
    {
        std::cerr << "no fmm configuration provided" << std::endl;
        return -1;
    }
   
    std::vector<Element*> src_elements, tgt_elements;
    unsigned int exp_terms, loc_terms, max_cell_elements;
    bool src_eq_tgt;
    read_fmm_config(argv[2], exp_terms, loc_terms, max_cell_elements);
    read_point2d_elements(argv[1],src_elements,tgt_elements,src_eq_tgt);
    FMM2D fmm(src_elements,
              tgt_elements,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    
    fmm.calculate();
    for(std::vector<Element*>::const_iterator it = tgt_elements.begin(); 
            it !=tgt_elements.end(); 
            ++it)
    {
        std::cout << (*it)->get_target_value() << std::endl;
    }
    std::cout << "Hello World!" << std::endl;
    return 0;
}

