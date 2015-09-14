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
#include "fmm2d_ada.h"
#include "gmres.h"
#include <cstdlib>
#include <set>
#include <iomanip>

#define POINT 0
#define POINT_ADA 0
#define CONST_EL 1
#define CONST_EL_ADA 0
#define OUTPUT_FMM 0
#define OUTPUT_COMP 0
#define FMM_GMRES 1
#define DIRECT_GMRES 1
#define DIRECT_SOLVE 0

#define TIMING 1


/*
 * For getRealTime()
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */
#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>	/* POSIX flags */
#include <time.h>	/* clock_gettime(), time() */
#include <sys/time.h>	/* gethrtime(), gettimeofday() */

#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#else
#error "Unable to define getRealTime( ) for an unknown OS."
#endif

/**
 * Returns the real time, in seconds, or -1.0 if an error occurred.
 *
 * Time is measured since an arbitrary and OS-dependent start time.
 * The returned real time is only useful for computing an elapsed time
 * between two calls to this function.
 */
double getRealTime( )
{
#if defined(_WIN32)
	FILETIME tm;
	ULONGLONG t;
#if defined(NTDDI_WIN8) && NTDDI_VERSION >= NTDDI_WIN8
	/* Windows 8, Windows Server 2012 and later. ---------------- */
	GetSystemTimePreciseAsFileTime( &tm );
#else
	/* Windows 2000 and later. ---------------------------------- */
	GetSystemTimeAsFileTime( &tm );
#endif
	t = ((ULONGLONG)tm.dwHighDateTime << 32) | (ULONGLONG)tm.dwLowDateTime;
	return (double)t / 10000000.0;

#elif (defined(__hpux) || defined(hpux)) || ((defined(__sun__) || defined(__sun) || defined(sun)) && (defined(__SVR4) || defined(__svr4__)))
	/* HP-UX, Solaris. ------------------------------------------ */
	return (double)gethrtime( ) / 1000000000.0;

#elif defined(__MACH__) && defined(__APPLE__)
	/* OSX. ----------------------------------------------------- */
	static double timeConvert = 0.0;
	if ( timeConvert == 0.0 )
	{
		mach_timebase_info_data_t timeBase;
		(void)mach_timebase_info( &timeBase );
		timeConvert = (double)timeBase.numer /
			(double)timeBase.denom /
			1000000000.0;
	}
	return (double)mach_absolute_time( ) * timeConvert;

#elif defined(_POSIX_VERSION)
	/* POSIX. --------------------------------------------------- */
#if defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
	{
		struct timespec ts;
#if defined(CLOCK_MONOTONIC_PRECISE)
		/* BSD. --------------------------------------------- */
		const clockid_t id = CLOCK_MONOTONIC_PRECISE;
#elif defined(CLOCK_MONOTONIC_RAW)
		/* Linux. ------------------------------------------- */
		const clockid_t id = CLOCK_MONOTONIC_RAW;
#elif defined(CLOCK_HIGHRES)
		/* Solaris. ----------------------------------------- */
		const clockid_t id = CLOCK_HIGHRES;
#elif defined(CLOCK_MONOTONIC)
		/* AIX, BSD, Linux, POSIX, Solaris. ----------------- */
		const clockid_t id = CLOCK_MONOTONIC;
#elif defined(CLOCK_REALTIME)
		/* AIX, BSD, HP-UX, Linux, POSIX. ------------------- */
		const clockid_t id = CLOCK_REALTIME;
#else
		const clockid_t id = (clockid_t)-1;	/* Unknown. */
#endif /* CLOCK_* */
		if ( id != (clockid_t)-1 && clock_gettime( id, &ts ) != -1 )
			return (double)ts.tv_sec +
				(double)ts.tv_nsec / 1000000000.0;
		/* Fall thru. */
	}
#endif /* _POSIX_TIMERS */

	/* AIX, BSD, Cygwin, HP-UX, Linux, OSX, POSIX, Solaris. ----- */
	struct timeval tm;
	gettimeofday( &tm, NULL );
	return (double)tm.tv_sec + (double)tm.tv_usec / 1000000.0;
#else
	return -1.0;		/* Failed. */
#endif
}

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

std::vector<double> direct_method_points(char* element_data, 
                                         std::vector<unsigned int> tgt_idx = std::vector<unsigned int>())
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
    
    //determine whether to compute everything or chosen target indexes
    if (tgt_idx.empty())
    {
        tgt_val.resize(tgt_x.size(),0);
    }
    else
    {
        tgt_val.resize(tgt_idx.size(),0);
    }
    
    

#if TIMING
	double time_start = getRealTime();
#endif
    //compute all interactions directly
    unsigned int num_tgt;
    unsigned int num_src = src_val.size();
    
    if(tgt_idx.empty())
    {
        num_tgt = tgt_x.size();
        for (unsigned int tgt = 0; tgt < num_tgt; tgt++)
        {
            double val = 0;
            double tgtx = tgt_x[tgt];
            double tgty = tgt_y[tgt];
            double diffx = 0;
            double diffy = 0;
            double dist = 0;
            for (unsigned int src = 0; src < num_src; src++)
            {
                diffx = tgtx - src_x[src];
                diffy = tgty - src_y[src];
                dist = std::sqrt(diffx * diffx + diffy * diffy);
                if (dist > 0)
                {
                    val += -src_val[src] * std::log(dist);
                }
            }
            tgt_val[tgt] = val;
        }
    }
    else //only compute at chosen indexes
    {
        num_tgt = tgt_idx.size();
        for (unsigned int idx = 0; idx < num_tgt; idx++)
        {
            double val = 0;
            double tgtx = tgt_x[tgt_idx[idx]];
            double tgty = tgt_y[tgt_idx[idx]];
            double diffx = 0;
            double diffy = 0;
            double dist = 0;
            for (unsigned int src = 0; src < num_src; src++)
            {
                diffx = tgtx - src_x[src];
                diffy = tgty - src_y[src];
                dist = std::sqrt(diffx * diffx + diffy * diffy);
                if (dist > 0)
                {
                    val += -src_val[src] * std::log(dist);
                }
            }
            tgt_val[idx] = val;
        }
    }
#if TIMING
    double time_end = getRealTime();
    double time = time_end-time_start;
    if(tgt_idx.size()) 
    {
        //scale up
        time *= tgt_x.size()/tgt_idx.size();
    }
    
    std::cout << "Direct method with " << num_src << " sources and " 
            << tgt_x.size() << " targets took " << time << " seconds "
			<< std::endl;
    
#endif
    
    return tgt_val;
}

void direct_method_elements(std::vector<Element*> const & src_el,
                            std::vector<Element*> const & tgt_el,
                            std::vector<unsigned int> tgt_idx = std::vector<unsigned int>())
{
    unsigned int num_src_el = src_el.size();
    unsigned int num_tgt_el;
    KernLapConstEl2D kernel;
    
#if TIMING
    double time_start = getRealTime();
#endif
    if(tgt_idx.empty())
    {
        num_tgt_el = tgt_el.size();
        for (unsigned int tgt = 0; tgt < num_tgt_el; tgt++)
        {
            double val = 0;
            Element* t = tgt_el[tgt];
                for (unsigned int src = 0; src < num_src_el; src++)
                {
                    Element* s = src_el[src];
                val += s->get_value() * kernel.direct(*t, *s);
                }
            t->set_target_value(val);
        }
    }
    else
    {
        num_tgt_el = tgt_idx.size();
        for(unsigned int idx = 0; idx < num_tgt_el; idx++)
        {
            double val = 0;
            Element* t = tgt_el[tgt_idx[idx]];
            for (unsigned int src = 0; src < num_src_el; src++)
            {
                Element* s = src_el[src];
                val += s->get_value() * kernel.direct(*t, *s);
            }
            t->set_target_value(val);
        }
    }
    
#if TIMING
    double time_end = getRealTime();
    double time = time_end-time_start;
    if(tgt_idx.size())
    {
        //scale up
        time*= tgt_el.size()/tgt_idx.size();
    }
    
    std::cout << "Direct method with " << num_src_el << " sources and " 
            << tgt_el.size() << " targets took " << time << " seconds"
            << std::endl;
#endif
}

void direct_method_el_coeff_mat(std::vector<Element*> const & src_el,
                                std::vector<Element*> const & tgt_el,
                                arma::mat & A_mat)
{
    assert(src_el.size() == tgt_el.size());
    unsigned int num_el = src_el.size();
    A_mat.resize(num_el,num_el);
    KernLapConstEl2D kernel;
    for(unsigned int i =0; i<num_el; i++)
    {
        Element *t = tgt_el[i];
        for(unsigned int j = 0; j<num_el; j++)
        {
            A_mat(i,j) = kernel.direct(*t,*(src_el[j]));
        }
    }
}

void point_fmm(std::vector<Element*> const & src_el,
               std::vector<Element*> const & tgt_el,
               int exp_terms,
               int loc_terms,
               int max_cell_elements,
               bool src_eq_tgt)
{
    FMM2D fmm(src_el,
              tgt_el,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    Laplace2DKernel laplace;
    fmm.set_kernel(laplace);
    
#if TIMING
    double time_start = getRealTime();
#endif
    fmm.calculate();
#if TIMING 
    double time_end = getRealTime();
    
    std::cout << "FMM with " << src_el.size() << " sources and " 
            << tgt_el.size() << " targets took " << time_end - time_start << " seconds "
			<< std::endl;
#endif
}

void point_ada_fmm(std::vector<Element*> const & src_el,
               std::vector<Element*> const & tgt_el,
               int exp_terms,
               int loc_terms,
               int max_cell_elements,
               bool src_eq_tgt)
{
    FMM2D_ADA fmm_ada(src_el,
              tgt_el,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    Laplace2DKernel laplace_ada;
    fmm_ada.set_kernel(laplace_ada);
    
#if TIMING
    double time_start_ada = getRealTime();
#endif
    fmm_ada.calculate();
#if TIMING
    double time_end_ada = getRealTime();
 std::cout << "Adaptive FMM with " << src_el.size() << " sources and " 
            << tgt_el.size() << " targets took " << time_end_ada - time_start_ada << " seconds "
			<< std::endl;
#endif
}

void const_el_fmm(std::vector<Element*> const & src_el,
               std::vector<Element*> const & tgt_el,
               int exp_terms,
               int loc_terms,
               int max_cell_elements,
               bool src_eq_tgt)
{
    FMM2D fmm(src_el,
              src_el,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    KernLapConstEl2D const_lap2d;
    fmm.set_kernel(const_lap2d);
#if TIMING
    double time_start = getRealTime();
#endif
    fmm.calculate();
#if TIMING 
    double time_end = getRealTime();
    
    std::cout << "FMM with " << src_el.size() << " sources and " 
            << tgt_el.size() << " targets took " << time_end - time_start << " seconds "
			<< std::endl;
#endif
}

void const_el_ada_fmm(std::vector<Element*> const & src_el,
               std::vector<Element*> const & tgt_el,
               int exp_terms,
               int loc_terms,
               int max_cell_elements,
               bool src_eq_tgt)
{
    FMM2D_ADA fmm_ada(src_el,
              tgt_el,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
    KernLapConstEl2D const_lap2d;
    fmm_ada.set_kernel(const_lap2d);
#if TIMING
    double time_start_ada = getRealTime();
#endif
    fmm_ada.calculate();
#if TIMING
    double time_end_ada = getRealTime();
    std::cout << "Adaptive FMM with " << src_el.size() << " sources and "
            << tgt_el.size() << " targets took " << time_end_ada - time_start_ada << " seconds "
            << std::endl;
#endif
}

std::vector<unsigned int> random_idx_choice(unsigned int range, unsigned int num)
{
    assert(range > 10*num);
    std::set<unsigned int> chosen_idx;
    double time = getRealTime();
    time = (time > 0) ? time : -time;
    unsigned int seed = round(time);
    srand(seed);
    while (chosen_idx.size() < num)
    {
        unsigned int choice = rand() % range;
        chosen_idx.insert(choice);
    }
    return std::vector<unsigned int>(chosen_idx.begin(),chosen_idx.end());
}

double mean_square_error(std::vector<double> const & correct_data, 
                         std::vector<double> const & cmp_data)
{
    assert(correct_data.size() == cmp_data.size());
    unsigned int data_size = correct_data.size();
    double sum_squares = 0;
    double sum_squares_diff = 0;
    for(unsigned int i = 0; i<correct_data.size(); i++)
    {
        double diff = correct_data[i] - cmp_data[i];
        diff *= diff;
        sum_squares_diff += diff;
        sum_squares += correct_data[i]*correct_data[i];
    }
    return std::sqrt(sum_squares_diff)/std::sqrt(sum_squares);
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
    
    std::vector<double> fmm1_res;
    std::vector<double> fmm_ada_res;
    
#if POINT || POINT_ADA
    read_point2d_elements(argv[1],src_elements,tgt_elements,src_eq_tgt);
#endif
    
#if POINT
    point_fmm(src_elements,
              tgt_elements,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
#if OUTPUT_COMP
    fmm1_res.resize(tgt_elements.size(),0);
    for(int i = 0; i<tgt_elements.size(); i++)
    {
        fmm1_res[i] = tgt_elements[i]->get_target_value();
    }
#endif
#endif

#if POINT_ADA
    point_ada_fmm(src_elements,
                  tgt_elements,
                  exp_terms,
                  loc_terms,
                  max_cell_elements,
                  src_eq_tgt);
#if OUTPUT_COMP
    fmm_ada_res.resize(tgt_elements.size(),0);
    for(int i = 0; i<tgt_elements.size(); i++)
    {
        fmm_ada_res[i] = tgt_elements[i]->get_target_value();
    }
#endif
#endif
    
#if POINT || POINT_ADA
    
    std::vector<double> direct_val;
    std::vector<unsigned int> tgt_idxes;
    if(src_elements.size() > 10000 && tgt_elements.size() > 1000)
    {
        // don't compute at all targets, make a random choice of 100 targets
        tgt_idxes = random_idx_choice(tgt_elements.size(),100);
        direct_val = direct_method_points(argv[1],tgt_idxes);
    }
    else
    {
        direct_val = direct_method_points(argv[1]);
    }
#endif
    
#if CONST_EL || CONST_EL_ADA
    read_const2d_elements(argv[1],src_elements,tgt_elements,src_eq_tgt);
#endif
    
#if CONST_EL
    const_el_fmm(src_elements,
                 tgt_elements,
                 exp_terms,
                 loc_terms,
                 max_cell_elements,
                 src_eq_tgt);
#if OUTPUT_COMP
    fmm1_res.resize(tgt_elements.size(),0);
    for(int i = 0; i<tgt_elements.size(); i++)
    {
        fmm1_res[i] = tgt_elements[i]->get_target_value();
    }
#endif
#endif
    
#if CONST_EL_ADA
    const_el_ada_fmm(src_elements,
                     tgt_elements,
                     exp_terms,
                     loc_terms,
                     max_cell_elements,
                     src_eq_tgt);
#if OUTPUT_COMP
    fmm_ada_res.resize(tgt_elements.size(),0);
    for(int i = 0; i<tgt_elements.size(); i++)
    {
        fmm_ada_res[i] = tgt_elements[i]->get_target_value();
    }
#endif
#endif

#if OUTPUT_FMM

    for(std::vector<Element*>::const_iterator it = tgt_elements.begin(); 
            it !=tgt_elements.end(); 
            ++it)
    {
        std::cout << (*it)->get_target_value() << std::endl;
    }
#endif
    
    
#if CONST_EL || CONST_EL_ADA
    std::vector<double> direct_val;
    std::vector<unsigned int> tgt_idxes;
    if(src_elements.size() > 10000 && tgt_elements.size() > 1000)
    {
        // don't compute at all targets, make a random choice of 100 targets
        tgt_idxes = random_idx_choice(tgt_elements.size(),100);
        direct_method_elements(src_elements,tgt_elements,tgt_idxes);
        //copy results
        direct_val.resize(tgt_idxes.size(),0);
        for(unsigned int i = 0; i<tgt_idxes.size(); i++)
        {
            direct_val[i] = tgt_elements[tgt_idxes[i]]->get_target_value();
        }
    }
    else
    {
        direct_method_elements(src_elements,tgt_elements);
        //copy results
        direct_val.resize(tgt_elements.size(),0);
        for(unsigned int i = 0; i < tgt_elements.size(); i++)
        {
            direct_val[i] = tgt_elements[i]->get_target_value();
        }
    }
#endif
    

#if (CONST_EL || CONST_EL_ADA) && (!CONST_EL || !CONST_EL_ADA) && FMM_GMRES
    std::vector<double> b_goals(tgt_elements.size(),1);
    std::vector<double> init_guess(tgt_elements.size(),0.5);
    std::vector<double> fmm_gmres_solution(tgt_elements.size(),0);
#if CONST_EL_ADA
    FMM2D_ADA fmm(src_elements,
              tgt_elements,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
#elif CONST_EL
    FMM2D fmm(src_elements,
              tgt_elements,
              exp_terms,
              loc_terms,
              max_cell_elements,
              src_eq_tgt);
#endif
    
    KernLapConstEl2D const_lap2d;
    fmm.set_kernel(const_lap2d);;

    FMM_GMRES_Solver fmm_solv(fmm,b_goals,init_guess,fmm_gmres_solution);
    double tol = 1e-8;
    int max_iter = 100;
    double fmm_gmres_start = getRealTime();
    fmm_solv.solve(max_iter,15,tol);
    double fmm_gmres_end = getRealTime();
    
    std::cout << "FMMGRES found solution after " << fmm_gmres_end-fmm_gmres_start 
            << " seconds and " << max_iter << " iterations "
            << "with tolerance " << tol << std::endl;
    std::cout <<"===" << std::endl;
    for(int i = 0; i<fmm_gmres_solution.size(); i++)
    {
        std::cout << fmm_gmres_solution[i] << std::endl;
    }
    std::cout<<"===" << std::endl;
    for(int i = 0; i<fmm_gmres_solution.size(); i++)
    {
        src_elements[i]->set_value(fmm_gmres_solution[i]);
    }
    fmm.recalculate();
    std::cout << "after recalc" << std::endl;
    for(int i = 0; i<fmm_gmres_solution.size(); i++)
    {
        std::cout << tgt_elements[i]->get_target_value() << std::endl;
    }
    std::cout << "===" << std::endl;
    
#endif
    
#if (DIRECT_SOLVE || DIRECT_GMRES) && (CONST_EL || CONST_EL_ADA) && (!CONST_EL || !CONST_EL_ADA) && FMM_GMRES
    arma::mat A_matrix;
    direct_method_el_coeff_mat(src_elements,tgt_elements,A_matrix);
#if DIRECT_SOLVE
    arma::vec b_dir(b_goals);
    arma::vec sol_dir;
    //std::cout << "Armadillo bug prevention: " << sol_dir << std::endl;
    double time_dir_solve_start = getRealTime();
    arma::solve(sol_dir,A_matrix,b_dir);
    double time_dir_solve_end = getRealTime();
    std::cout << "direct solution after " <<time_dir_solve_end - time_dir_solve_start 
            << " seconds"<< std::endl;
    std::cout << "===" << std::endl;
    std::cout << sol_dir << std::endl;
    std::cout << "===" << std::endl;
#endif
    
#if DIRECT_GMRES
    int restart_m = 15;
    int max_iter_dir = 100;
    double tol_dir = 1e-8;
    arma::vec b_gmres(b_goals);
    arma::vec x;
    x = arma::conv_to<arma::vec>::from(init_guess);
    arma::mat H(restart_m+1,restart_m+1);
    NoPrecond no_pre;
    std::cout << "Armadillo bug prevention: " << x << std::endl;
    double time_dir_solve_gmres_start = getRealTime();
    GMRES<arma::mat,
          arma::vec,
          NoPrecond,
          arma::mat,
          double>(A_matrix,x,b_gmres,no_pre,H,restart_m,max_iter_dir,tol_dir);
    double time_dir_solve_gmres_end = getRealTime();
    std::cout << "got direct GMRES solution after " << time_dir_solve_gmres_end - time_dir_solve_gmres_start
            << " seconds and " << max_iter << " iterations "
            << "with tolerance " << tol << std::endl;
    std::cout << "===" << std::endl;
    std::cout << std::setprecision(15) << x << std::endl;
    std::cout << "===" << std::endl;
    std::cout << std::setprecision(15) << A_matrix*x << std::endl;
    std::cout << "===" << std::endl;
#endif
    
#if DIRECT_SOLVE
    double mean_square_gmres_dir = mean_square_error(
                                        arma::conv_to<std::vector<double> >::from(sol_dir),
                                        fmm_gmres_solution);
    std::cout << "Mean Square Error of direct solution vs FMMGMRES is " << mean_square_gmres_dir << std::endl;
#endif
    
#if DIRECT_GMRES
    double mean_square_gmres_gmres = mean_square_error(
                                        arma::conv_to<std::vector<double> >::from(x),
                                        fmm_gmres_solution);
    std::cout << "Mean Square Error of GMRES solution vs FMMGMRES is " << mean_square_gmres_gmres << std::endl;        
#endif
    
#endif
    
    //calculate errors and print them
    
#if OUTPUT_COMP
    /* std::cout << direct_val.size() << std::endl;
    unsigned int num_tgts = direct_val.size();
    for(int i = 0; i<num_tgts; i++)
    {
        std::cout << "diff: dir: " << direct_val[i] <<" fmm: " << tgt_elements[i]->get_target_value() << " diff: " << direct_val[i]-tgt_elements[i]->get_target_value() << std::endl;
    }
    */
    
    std::vector<double> res1;
    std::vector<double> res_ada;
    
    if(tgt_idxes.empty())
    {
       res1 = fmm1_res;
       res_ada = fmm_ada_res;
    } 
    else
    {
        if(!fmm1_res.empty())
        {
            for(int i = 0; i<tgt_idxes.size(); i++)
            {
                res1.push_back(fmm1_res[tgt_idxes[i]]);
            }
        }
        if(!fmm_ada_res.empty())
        {
            for(int i = 0; i<tgt_idxes.size(); i++)
            {
                res_ada.push_back(fmm_ada_res[tgt_idxes[i]]);
            }
        }
    }
    if (!fmm1_res.empty())
    {
        double msq_error_1 = mean_square_error(direct_val, res1);
        std::cout << "Mean Square Error of FMM vs. direct is " << std::setprecision(15)
                << msq_error_1 << std::endl;
    }
    if (!fmm_ada_res.empty())
    {
        double msq_error_ada = mean_square_error(direct_val, res_ada);
        std::cout << "Mean Square Error of adaptive FMM vs. direct is " << std::setprecision(15)
                << msq_error_ada << std::endl;
    }
    
#endif
    std::cout << "Hello World!" << std::endl;
    return 0;
}

