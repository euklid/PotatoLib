#include <iostream>
#include "tree.h"
#include "element.h"
#include "point_element.h"
#include "point.h"
#include "fmm2d.h"
#include "kernel_laplace_point_2d.h"
#include "fmm_gmres_solver.h"
#include "constant_element_2d.h"
#include "kernel_laplace_constant_element_2d.h"
#include <time.h>



#define POINT 0
#define CONST_EL 1
#define OUTPUT_FMM 1
#define OUTPUT_COMP 0
#define FMM_GMRES 1
#define DIRECT_GMRES 0
#define DIRECT_SOLVE 0

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

//FIXME: be cool and use JSON for data configuration

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

void read_const2d_elements(char* filename,
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
     * #src_eq_tgt(0|1) 
     * 1
     * # num nodes 
     * 4
     * # node_nr_ascending(int)  node_x(double) node_y(double)
     * 0 0.0 0.0
     * 1 0.0 1.0
     * 2 1.0 1.0
     * 3 0.0 1.0
     * 
     * #start and end nodes of Elements. We assume that the enumeration of
     * #elements is counter-clock-wise
     * # node_start_nr node_end_nr s(0|1) t(0|1) init_val
     * 0 1 1 1 0.4
     * 1 2 1 1 0.4
     * 2 3 1 1 0.4
     * 3 0 1 1 0.4
     * 
     */
    
    file >> src_eq_tgt;
    double x, y, init_val;
    int source, target;
    int id = 0;
    int num_nodes;
    file >> num_nodes;
    std::vector<Point> id_nodes(num_nodes,Point(2));
    for(int i = 0; i< num_nodes; i++)
    {
        file >> id >> x >> y;
        id_nodes[i][0] = x;
        id_nodes[i][1] = y;
    }
    id = 0;
    int start_id, end_id;
    while(file >> start_id >> end_id >> source >> target >> init_val)
    {
        int type = source*Element::SOURCE + target*Element::TARGET;
        ConstEl2D* el = new ConstEl2D(id_nodes[start_id],id_nodes[end_id],id,type);
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

std::vector<double> direct_method_points(char* element_data)
{
    //read data
    std::ifstream file(element_data);
    std::vector<double> src_x, src_y, src_val;
    std::vector<double> tgt_x, tgt_y, tgt_val;
    src_x.reserve(1000000);
    src_y.reserve(1000000);
    src_val.reserve(1000000);
    tgt_x.reserve(1000000);
    tgt_y.reserve(1000000);// We've got PLENTY of RAM!
    double src_eq_tgt;
    assert(file.is_open());
    file >> src_eq_tgt;
    double x,y,init_val;
    int source, target;
    while(file >> x >> y >> source >> target >> init_val)
    {
        if(source)
        {
            src_x.push_back(x);
            src_y.push_back(y);
            src_val.push_back(init_val);
        }
        if(target)
        {
            tgt_x.push_back(x);
            tgt_y.push_back(y);
        }
    }
    
    file.close();
    tgt_val.resize(tgt_x.size(),0);
    
    timespec  time_dir_s, time_dir_e;
    time_t  dir_sec;
    long  dir_n_sec;
  
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_dir_s);
    
    //compute all interactions directly
    unsigned int num_tgt = tgt_x.size();
    unsigned int num_src = src_val.size();
    for(unsigned int tgt = 0; tgt<num_tgt; tgt++)
    {
        double val = 0;
        double tgtx = tgt_x[tgt];
        double tgty = tgt_y[tgt];
        double diffx = 0;
        double diffy = 0;
        double dist = 0;
        for(unsigned int src = 0; src<num_src; src++)
        {
            diffx = tgtx-src_x[src];
            diffy = tgty-src_y[src];
            dist = std::sqrt(diffx*diffx+diffy*diffy);
            if(dist > 0)
            {
                val += -src_val[src]*std::log(dist);
            }
        }
        tgt_val[tgt] = val;
    }
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_dir_e);
    dir_sec = time_dir_e.tv_sec-time_dir_s.tv_sec;
    dir_n_sec = std::max(time_dir_e.tv_nsec - time_dir_s.tv_nsec,
                         1000000000l - (time_dir_e.tv_nsec - time_dir_s.tv_nsec)) ;
    
    std::cout << "Direct method with " << num_src << " sources and " 
            << num_tgt << " targets took " << dir_sec << " seconds and "
            << dir_n_sec << " nanoseconds" << std::endl;
    
    
    return tgt_val;
}

void direct_method_elements(std::vector<Element*> const & src_el,
                            std::vector<Element*> const & tgt_el)
{
    unsigned int num_src_el = src_el.size();
    unsigned int num_tgt_el = tgt_el.size();
    KernLapConstEl2D kernel;
    
    timespec  time_dir_s, time_dir_e;
    time_t  dir_sec;
    long  dir_n_sec;
  
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_dir_s);
    
    for(unsigned int tgt = 0; tgt < num_tgt_el; tgt++)
    {
        double val = 0;
        Element* t = tgt_el[tgt];
        for(unsigned int src = 0; src < num_src_el; src++)
        {
            Element* s = src_el[src];
            val += s->get_value()*kernel.direct(*t,*s);
        }
        t->set_target_value(val);
    }
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_dir_e);
    dir_sec = time_dir_e.tv_sec-time_dir_s.tv_sec;
    dir_n_sec = std::max(time_dir_e.tv_nsec - time_dir_s.tv_nsec,
                         1000000000l - (time_dir_e.tv_nsec - time_dir_s.tv_nsec)) ;
    
    std::cout << "Direct method with " << num_src_el << " sources and " 
            << num_tgt_el << " targets took " << dir_sec << " seconds and "
            << dir_n_sec << " nanoseconds" << std::endl;
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
   
    struct timespec time_fmm_s, time_fmm_e, time_dir_s, time_dir_e;
    time_t fmm_sec, dir_sec;
    long fmm_n_sec, dir_n_sec;
    
    std::vector<Element*> src_elements, tgt_elements;
    unsigned int exp_terms, loc_terms, max_cell_elements;
    bool src_eq_tgt;
    read_fmm_config(argv[2], exp_terms, loc_terms, max_cell_elements);
            
#if POINT
    read_point2d_elements(argv[1],src_elements,tgt_elements,src_eq_tgt);
    
    
    FMM2D fmm(src_elements,
              tgt_elements,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    Laplace2DKernel laplace;
    fmm.set_kernel(laplace);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_fmm_s);
    fmm.calculate();
#endif
    
#if CONST_EL
    read_const2d_elements(argv[1],src_elements,tgt_elements,src_eq_tgt);
    FMM2D fmm(src_elements,
              tgt_elements,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    KernLapConstEl2D const_lap2d;
    fmm.set_kernel(const_lap2d);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_fmm_s);
    fmm.calculate();
    
#endif
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_fmm_e);
    fmm_sec = time_fmm_e.tv_sec-time_fmm_s.tv_sec;
    fmm_n_sec = std::max(time_fmm_e.tv_nsec - time_fmm_s.tv_nsec,1000000000l - time_fmm_e.tv_nsec - time_fmm_s.tv_nsec);
    
    std::cout << "FMM with " << src_elements.size() << " sources and " 
            << tgt_elements.size() << " targets took " << fmm_sec << " seconds and "
            << fmm_n_sec << " nanoseconds" << std::endl;
    

#if OUTPUT_FMM

    for(std::vector<Element*>::const_iterator it = tgt_elements.begin(); 
            it !=tgt_elements.end(); 
            ++it)
    {
        std::cout << (*it)->get_target_value() << std::endl;
    }
    fmm.recalculate();
    for(std::vector<Element*>::const_iterator it = tgt_elements.begin(); 
            it !=tgt_elements.end(); 
            ++it)
    {
        std::cout << (*it)->get_target_value() << std::endl;
    }
#endif
    
std::vector<double> direct_val;

#if POINT
    // direct method to compare
    direct_val = direct_method_points(argv[1]);
#endif
    
#if CONST_EL
    unsigned int num_tgt_el = tgt_elements.size();
    direct_val.resize(num_tgt_el,0);
    for(unsigned int i = 0; i<num_tgt_el; i++)
    {
        direct_val[i] = tgt_elements[i]->get_target_value();
    }
    direct_method_elements(src_elements,tgt_elements);

#endif
    
    //calculate errors and print them

#if CONST_EL && FMM_GMRES
    std::vector<double> b_goals(tgt_elements.size(),1);
    std::vector<double> init_guess(tgt_elements.size(),0.0001);
    std::vector<double> solution(tgt_elements.size(),0);
    FMM_GMRES_Solver fmm_solv(fmm,b_goals,init_guess,solution);
    double tol = 1e-8;
    fmm_solv.solve(100,100,tol);
    
    for(int i = 0; i<solution.size(); i++)
    {
        std::cout << solution[i] << std::endl;
    }
    for(int i = 0; i<solution.size(); i++)
    {
        src_elements[i]->set_value(solution[i]);
    }
    fmm.recalculate();
    std::cout << "after recalc" << std::endl;
    for(int i = 0; i<solution.size(); i++)
    {
        std::cout << tgt_elements[i]->get_target_value() << std::endl;
    }
    
#endif
    
#if OUTPUT_COMP
    std::cout << direct_val.size() << std::endl;
    unsigned int num_tgts = direct_val.size();
    for(int i = 0; i<num_tgts; i++)
    {
        std::cout << "diff: dir: " << direct_val[i] <<" fmm: " << tgt_elements[i]->get_target_value() << " diff: " << direct_val[i]-tgt_elements[i]->get_target_value() << std::endl;
    }
#endif
    std::cout << "Hello World!" << std::endl;
    return 0;
}

