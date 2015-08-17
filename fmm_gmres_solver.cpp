#include "fmm_gmres_solver.h"
#include "gmres.h"
#include <armadillo>

class Operator
{
public:
    Operator(FMM & fmm) :
    m_fmm(fmm)
    {}
    
    arma::vec operator*(arma::vec vec) const
    {
        std::vector<Element*> const & src_el = m_fmm.get_src_elements();
        std::vector<Element*> const & tgt_el = m_fmm.get_tgt_elements();
        assert(vec.size() == src_el.size());
        assert(tgt_el.size() == src_el.size());
        unsigned int num_el = src_el.size();
        unsigned int num_tgt_el = tgt_el.size();
        for (unsigned int i = 0; i<num_el; i++) {
            src_el[i]->set_value(vec[i]);
        }
        m_fmm.recalculate();
        arma::vec res(num_tgt_el);
        for (unsigned int i = 0; i<num_tgt_el; i++) {
            res[i] = tgt_el[i]->get_target_value();
        }
        return res;
    }
    
private:
    FMM & m_fmm;
};

class ArmaMatWrapper : public arma::mat
{
public:
    ArmaMatWrapper(arma::mat & mat)
    : arma::mat(mat)
    {}
    
    arma::vec solve(arma::vec const & vec) const
    {
        arma::vec sol = arma::solve(static_cast<arma::mat>(*this), vec);
        return sol;
    }
};

class ArmaSpMatWrapper : public arma::sp_mat
{
public:
    ArmaSpMatWrapper(arma::sp_mat & mat)
    : arma::sp_mat(mat)
    {}
    
    arma::vec solve(arma::vec const & vec) const
    {
        arma::vec sol = arma::spsolve(static_cast<arma::sp_mat>(*this), vec);
        return sol;
    }
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
                                   std::vector<double> const & boundary_goals,
                                   std::vector<double> const & solution) :
    m_fmm(fmm),
    m_boundary_goals(boundary_goals),
    m_solution(solution)
{}

void FMM_GMRES_Solver::solve(int max_iterations, int m, double &tolerance)
{
    Operator A(m_fmm);
    arma::vec x=arma::ones(m_solution.size());
    arma::vec b(m_boundary_goals);
    arma::mat _H;
    arma::sp_mat _M;
    ArmaMatWrapper H(_H);
    ArmaSpMatWrapper M(_M);
    
    GMRES<Operator,
          arma::vec,
          ArmaSpMatWrapper,
          ArmaMatWrapper,
          double>(A, x, b, M, H, m, max_iterations, tolerance);
    
    
    
}
