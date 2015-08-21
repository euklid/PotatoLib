#include "fmm_gmres_solver.h"
#include "gmres.h"
#include <armadillo>

class Operator
{
public:
    Operator(FMM & fmm) :
    m_fmm(fmm), m_has_precond(false)
    {}
    
    arma::vec operator*(arma::vec const & vec)
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

class Precond
{
public:
    Precond(FMM & fmm)
    : m_fmm(fmm),m_has_prec_data(false)
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
            m_prec = &m_fmm.get_precond();
            m_has_prec_data = true;
            block_LUP();
        }
        std::vector<unsigned int> const & permutation = m_fmm.get_element_permutation();
        assert(permutation.size() == vec.size());
        arma::vec perm_vec(vec);
        unsigned int perm_size = permutation.size();
        for(unsigned int i = 0; i< perm_size; i++)
        {
            perm_vec(i) = vec(permutation[i]); // i-th entry is permutation(i)-th element
        }
        arma::vec perm_sol;
        LUP_solve(perm_vec, perm_sol);
        arma::vec sol;
        
        //bring solution in correct order
        for(unsigned int i = 0; i<perm_size; i++)
        {
            sol(permutation[i]) = perm_sol(i);
        }
        return sol;
    }
    
private:
    void block_LUP()
    {
        unsigned int num_blocks = m_prec->size();
        m_L.reserve(num_blocks);
        m_U.reserve(num_blocks);
        m_P.resize(num_blocks,0);
        for(unsigned int i = 0; i<num_blocks; i++)
        {
            arma::mat L,U,P;
            arma::lu(m_L[i],m_U[i],P,m_prec->at(i));
            unsigned int block_size = P.n_cols;
            arma::uvec perm(block_size);
            for(unsigned int j = 0; j < block_size; j++)
            {
                perm(j) = j;
            }
            m_P[i] = arma::conv_to<arma::uvec >::from(P*perm);
        }
    }
    
    void LUP_solve(arma::vec const & vec, arma::vec & sol)
    {
        std::vector<unsigned int> const & block_starts = m_fmm.get_prec_block_starts();
        unsigned int num_blocks = block_starts.size();
        for(int i = 0; i<num_blocks; i++)
        {
            // Ax=b <==> LUx = Pb <==> Ly = Pb and y = Ux
            // therefore first solve for y and then solve for x, knowing L and U
            // are triangular matrices and armadillo can use this to accelerate
            // the solution
            
            unsigned int block_size = m_prec->at(i).size();
            arma::vec block_vec(block_size);
            unsigned int block_start = block_starts[i];
            assert(m_P[i].size() == block_size);
            
            //initialize Pb
            for(unsigned int j = 0; j<block_size; j++)
            {
                block_vec(j) = vec(block_start+m_P[i](j));
            }
            // solve for y and x and assign to solution
            arma::vec y,x;
            arma::auxlib::solve_tr(y,m_L[i],block_vec,1); //0: upper, 1:lower
            arma::auxlib::solve_tr(x,m_U[i],y,0);
            sol.subvec(block_start,block_start+block_size-1) = x;
        }
    }
    
private:
    FMM & m_fmm;
    std::vector<arma::mat> const * m_prec;
    bool m_has_prec_data;
    
    // block matrices composed of L resp. U matrices
    std::vector<arma::mat > m_L,m_U; 
    
    //permutations for each block for LU decomp.
    std::vector<arma::uvec > m_P; 
};

FMM_GMRES_Solver::FMM_GMRES_Solver(FMM & fmm,
                                   std::vector<double> & boundary_goals,
                                   std::vector<double> & init_guess,
                                   std::vector<double> & solution) :
    m_fmm(fmm),
    m_boundary_goals(boundary_goals),
    m_init_guess(init_guess),
    m_solution(solution)
{}

void FMM_GMRES_Solver::solve(int max_iterations, int m, double &tolerance)
{
    assert(m_fmm.get_src_elements().size() == m_fmm.get_tgt_elements().size());
    Operator A(m_fmm);
    arma::vec x=arma::conv_to<arma::vec>::from(m_init_guess);
    arma::vec b(m_boundary_goals);
    arma::mat H(m,m);
    Precond M(m_fmm);
    
    GMRES<Operator,
          arma::vec,
          Precond,
          arma::mat,
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