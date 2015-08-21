#include <iostream>
#include "tree.h"
#include "element.h"
#include "point_element.h"
#include "point.h"
#include "fmm2d.h"
#include "kernel_laplace_point_2d.h"
#include "fmm_gmres_solver.h"
#include <time.h>

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

std::vector<double> direct_method(char* element_data)
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
    Laplace2DKernel laplace;
    fmm.set_kernel(laplace);
    
    struct timespec time_fmm_s, time_fmm_e, time_dir_s, time_dir_e;
    time_t fmm_sec, dir_sec;
    long fmm_n_sec, dir_n_sec;
  
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_fmm_s);
    
    fmm.calculate();
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_fmm_e);
    fmm_sec = time_fmm_e.tv_sec-time_fmm_s.tv_sec;
    fmm_n_sec = std::max(time_fmm_e.tv_nsec - time_fmm_s.tv_nsec,1000000000l - time_fmm_e.tv_nsec - time_fmm_s.tv_nsec);
    
    std::cout << "FMM with " << src_elements.size() << " sources and " 
            << tgt_elements.size() << " targets took " << fmm_sec << " seconds and "
            << fmm_n_sec << " nanoseconds" << std::endl;
    
   /* 
    for(std::vector<Element*>::const_iterator it = tgt_elements.begin(); 
            it !=tgt_elements.end(); 
            ++it)
    {
        std::cout << (*it)->get_target_value() << std::endl;
    }
    */
    
    // direct method to compare
    //std::vector<double> direct_val = direct_method(argv[1]);
    
    //calculate errors and print them
    /*
    std::cout << direct_val.size() << std::endl;
    unsigned int num_tgts = direct_val.size();
    for(int i = 0; i<num_tgts; i++)
    {
        std::cout << "diff: dir: " << direct_val[i] <<" fmm: " << tgt_elements[i]->get_target_value() << std::endl;
    }
    */
    std::cout << "Hello World!" << std::endl;
    return 0;
}

