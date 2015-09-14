#ifndef FMM_GMRES_SOLVER_H
#define FMM_GMRES_SOLVER_H

#include "fmm.h"

class FMM_GMRES_Solver {
    
public:
    /**
     *  Use a FMM to calculate missing boundary values for given
     *  boundary conditions using GMRES
     *  Caveats: only use with source elements being target elements 
     *  simultaneously!
     */
    FMM_GMRES_Solver(FMM  & fmm,
                     std::vector<double> & boundary_goals,
                     std::vector<double> const & init_guess,
                     std::vector<double> & solution);
    
    /**
     *  solve with GMRES
     *  @param max_iterations maximal iterations of GMRES
     *  @param m using GMRES(m) with restart after m iterations
     *  @param tolerance for convergence, after execution the norm of 
     *  residual
     */
    void solve(int& max_iterations, int m, double & tolerance);
    
private:
    FMM & m_fmm;
    std::vector<double> & m_boundary_goals;
    std::vector<double> & m_solution;
    std::vector<double> const & m_init_guess;
};

class NoPrecond
{
public:
    NoPrecond(){}
    arma::vec solve(arma::vec const & vec)
    {
        return vec;
    }
};



#endif