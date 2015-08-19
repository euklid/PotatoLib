#include "fmm_gmres_solver.h"
#include "gmres.h"
#include <armadillo>

class Operator
{
public:
    Operator(FMM & fmm) :
    m_fmm(fmm), m_has_precond(false)
    {}
    
    arma::vec operator*(arma::vec vec)
    {
        assert(vec.size() == m_fmm.get_src_elements().size());
        std::vector<Element*> const & source_el = m_fmm.get_src_elements();
        unsigned int num_el = source_el.size();
        for (unsigned int i = 0; i<num_el; i++) 
        {
            source_el[i]->set_value(vec(i));
        }
        if(!m_has_precond)
        {
            m_has_precond = m_fmm.get_precond().size();
            if(!m_has_precond)
            {
                m_fmm.calculate(true);
            }
            else
            {
                m_fmm.recalculate();
            }
        } 
        else
        {
            m_fmm.recalculate();
        }
        arma::vec res(num_el);
        //ATTENTION: We assume that we have combined elements that are source
        //           and target at the same time
        for (unsigned int i = 0; i<num_el; i++) 
        {
            res(i) = source_el[i]->get_target_value();
        }
        return res;
    }
    
private:
    FMM & m_fmm;
    bool m_has_precond;
};

class ArmaMatWrapper : public arma::mat
{
public:
    ArmaMatWrapper(arma::mat & mat)
    : arma::mat(mat)
    {}
    
    ArmaMatWrapper(int rows, int cols)
    : arma::mat(rows,cols)
    {}
    
    arma::vec solve(arma::vec const & vec) const
    {
        arma::vec sol = arma::solve(static_cast<arma::mat>(*this), vec);
        return sol;
    }
};

class Precond : public arma::sp_mat
{
public:
    Precond(FMM & fmm)
    : arma::sp_mat(), m_fmm(fmm),m_has_prec_data(false)
    {}
    
    /**
     * solves system (this)x = vec and returns x
     * Further, because this is the preconditioner computed by a FMM,
     * we have to first solve the henn and egg problem to obtain the matrix
     * entries before we can actually do anything and this is only possible 
     * after the FMM calculation has been completed at least once.
     * KNOWING the preconditioner will be called AFTER the first calculation
     * of the FMM, get the preconditioning matrix from the FMM.
     * @param vec vector on right side
     * @return solution vector x
     */
    arma::vec solve(arma::vec const & vec)
    {
        if(!m_has_prec_data)
        {
            (dynamic_cast<arma::sp_mat*>(this))->operator =(m_fmm.get_precond());
            m_has_prec_data = true;
        }
        std::vector<unsigned int> const & permutation = m_fmm.get_element_permutation();
        assert(permutation.size() == vec.size());
        arma::vec perm_vec(vec);
        unsigned int perm_size = permutation.size();
        for(unsigned int i = 0; i< perm_size; i++)
        {
            perm_vec(i) = vec(permutation[i]); // i-th entry is permutation(i)-th element
        }
        arma::vec perm_sol = arma::spsolve(static_cast<arma::sp_mat>(*this), vec);
        arma::vec sol(perm_sol);
        
        //bring solution in correct order
        for(unsigned int i = 0; i<perm_size; i++)
        {
            sol(permutation[i]) = perm_sol(i);
        }
        return sol;
    }
private:
    FMM & m_fmm;
    bool m_has_prec_data;
};

//FIXME: how to switch fast between target/source values of elements and boundary conditions/goals.
// 1. idea: use a backing array for elements and override get_value() functions to access these
//          the backing array then would be used to solve the equations with armadillo
// 2. idea: be naive and copy
// 3. idea: use the arma types for internal vectors, what internal vectors?
//          * only used vector as C-array replacement, even have own Point class
//            + combine ideas 1 and 3 and use armadillo vectors as internal storage for values
//            + FMM still needs to produce a preconditioner

 

FMM_GMRES_Solver::FMM_GMRES_Solver(FMM & fmm,
                                   std::vector<double> & boundary_goals,
                                   std::vector<double> & solution) :
    m_fmm(fmm),
    m_boundary_goals(boundary_goals),
    m_solution(solution)
{}

void FMM_GMRES_Solver::solve(int max_iterations, int m, double &tolerance)
{
    Operator A(m_fmm);
    arma::vec x=arma::ones(m_solution.size());
    arma::vec b(m_boundary_goals);
    //arma::mat _H(m,m);
    //arma::sp_mat _M;
    ArmaMatWrapper H(m,m);
    Precond M(m_fmm);
    
    GMRES<Operator,
          arma::vec,
          Precond,
          ArmaMatWrapper,
          double>(A, x, b, M, H, m, max_iterations, tolerance);
    
    //output result
    for(unsigned int i = 0; i<m_boundary_goals.size(); i++)
    {
        m_boundary_goals[i] = b(i);
    }
    
    for(unsigned int i = 0; i<m_solution.size(); i++)
    {
        m_solution[i] = x(i);
    }
}
