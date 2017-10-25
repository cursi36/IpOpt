#ifndef __HS071_NLP_HPP__
#define __HS071_NLP_HPP__

#include <IpTNLP.hpp>
#include <cassert>
#include <iostream>

using namespace Ipopt;

class HS071_NLP: public TNLP
{
  
public:
  HS071_NLP()
  {}
  
  virtual ~HS071_NLP()
  {}
  
  
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  //n variables
  n = 4;
  
  //m constraints
  m = 2;
  
    // in this example the Jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 8;
  
    // the Hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 10;
  
    // use the C style indexing (0-based).the numbering style used for row/col entries in the sparse matrix format 
  index_style = TNLP::C_STYLE;

  return true;
}


// Give IPOPT the value of the bounds on the variables and constraints. 
/*
  
    n: (in), the number of variables in the problem (dimension of $ x$).
    x_l: (out) the lower bounds $ x^L$ for $ x$.
    x_u: (out) the upper bounds $ x^U$ for $ x$.
    m: (in), the number of constraints in the problem (dimension of $ g(x)$).
    g_l: (out) the lower bounds $ g^L$ for $ g(x)$.
    g_u: (out) the upper bounds $ g^U$ for $ g(x)$.
*/

virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u)
{
  n = 4;
  
  m = 2;
  
    for (Index i=0; i<4; i++)
    x_l[i] = 1.0;
    
  for (Index i=0; i<4; i++)
    x_u[i] = 5.0;
  
   g_l[0] = 25;
  
   g_l[1] = g_u[1] = 40.0;

  return true;
  
}


// returns the initial point for the problem
virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
 // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish to use a warmstart option
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
  
  
  // initialize to the given starting point
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;

  return true;
  
  
}

// returns the value of the objective function
virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 4);

  obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
//The gradient array is in the same order as the $ x$ variables (i.e., the gradient of the objective with respect to x[2] should be put in grad_f[2]). 
bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 4);

  grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
  grad_f[1] = x[0] * x[3];
  grad_f[2] = x[0] * x[3] + 1;
  grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

  return true;
}


// return the value of the constraints: g(x)
virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 4);
  assert(m == 2);

  g[0] = x[0] * x[1] * x[2] * x[3];
  g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];

  return true;
}

// return the structure or values of the jacobian

/*

    n: (in), the number of variables in the problem (dimension of $ x$).
    x: (in), the values for the primal variables, $ x$, at which the constraint Jacobian, $ \nabla g(x)^T$, is to be evaluated.
    new_x: (in), false if any evaluation method was previously called with the same values in x, true otherwise.
    m: (in), the number of constraints in the problem (dimension of $ g(x)$).
    n_ele_jac: (in), the number of nonzero elements in the Jacobian (dimension of iRow, jCol, and values).
    iRow: (out), the row indices of entries in the Jacobian of the constraints.
    jCol: (out), the column indices of entries in the Jacobian of the constraints.
    values: (out), the values of the entries in the Jacobian of the constraints.
*/

virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    // this particular jacobian is dense
    iRow[0] = 0;
    jCol[0] = 0;
    iRow[1] = 0;
    jCol[1] = 1;
    iRow[2] = 0;
    jCol[2] = 2;
    iRow[3] = 0;
    jCol[3] = 3;
    iRow[4] = 1;
    jCol[4] = 0;
    iRow[5] = 1;
    jCol[5] = 1;
    iRow[6] = 1;
    jCol[6] = 2;
    iRow[7] = 1;
    jCol[7] = 3;
  }
  else {
    // return the values of the jacobian of the constraints

    values[0] = x[1]*x[2]*x[3]; // 0,0
    values[1] = x[0]*x[2]*x[3]; // 0,1
    values[2] = x[0]*x[1]*x[3]; // 0,2
    values[3] = x[0]*x[1]*x[2]; // 0,3

    values[4] = 2*x[0]; // 1,0
    values[5] = 2*x[1]; // 1,1
    values[6] = 2*x[2]; // 1,2
    values[7] = 2*x[3]; // 1,3
  }

  return true;
}

//return the structure or values of the hessian
/*
 * 
    n: (in), the number of variables in the problem (dimension of $ x$).
    x: (in), the values for the primal variables, $ x$, at which the Hessian is to be evaluated.
    new_x: (in), false if any evaluation method was previously called with the same values in x, true otherwise.
    obj_factor: (in), factor in front of the objective term in the Hessian, $ \sigma_f$.
    m: (in), the number of constraints in the problem (dimension of $ g(x)$).
    lambda: (in), the values for the constraint multipliers, $ \lambda$, at which the Hessian is to be evaluated.
    new_lambda: (in), false if any evaluation method was previously called with the same values in lambda, true otherwise.
    nele_hess: (in), the number of nonzero elements in the Hessian (dimension of iRow, jCol, and values).
    iRow: (out), the row indices of entries in the Hessian.
    jCol: (out), the column indices of entries in the Hessian.
    values: (out), the values of the entries in the Hessian.
*/


virtual bool eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < 4; row++) {
      for (Index col = 0; col <= row; col++) {
        iRow[idx] = row;
        jCol[idx] = col;
        idx++;
      }
    }

    assert(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    // fill the objective portion
    values[0] = obj_factor * (2*x[3]); // 0,0

    values[1] = obj_factor * (x[3]);   // 1,0
    values[2] = 0.;                    // 1,1

    values[3] = obj_factor * (x[3]);   // 2,0
    values[4] = 0.;                    // 2,1
    values[5] = 0.;                    // 2,2

    values[6] = obj_factor * (2*x[0] + x[1] + x[2]); // 3,0
    values[7] = obj_factor * (x[0]);                 // 3,1
    values[8] = obj_factor * (x[0]);                 // 3,2
    values[9] = 0.;                                  // 3,3


    // add the portion for the first constraint
    values[1] += lambda[0] * (x[2] * x[3]); // 1,0

    values[3] += lambda[0] * (x[1] * x[3]); // 2,0
    values[4] += lambda[0] * (x[0] * x[3]); // 2,1

    values[6] += lambda[0] * (x[1] * x[2]); // 3,0
    values[7] += lambda[0] * (x[0] * x[2]); // 3,1
    values[8] += lambda[0] * (x[0] * x[1]); // 3,2

    // add the portion for the second constraint
    values[0] += lambda[1] * 2; // 0,0

    values[2] += lambda[1] * 2; // 1,1

    values[5] += lambda[1] * 2; // 2,2

    values[9] += lambda[1] * 2; // 3,3
  }

  return true;
}


/*
    status: (in), gives the status of the algorithm as specified in IpAlgTypes.hpp,
        SUCCESS: Algorithm terminated successfully at a locally optimal point, satisfying the convergence tolerances (can be specified by options).
        MAXITER_EXCEEDED: Maximum number of iterations exceeded (can be specified by an option).
        CPUTIME_EXCEEDED: Maximum number of CPU seconds exceeded (can be specified by an option).
        STOP_AT_TINY_STEP: Algorithm proceeds with very little progress.
        STOP_AT_ACCEPTABLE_POINT: Algorithm stopped at a point that was converged, not to ``desired'' tolerances, but to ``acceptable'' tolerances (see the acceptable-... options).
        LOCAL_INFEASIBILITY: Algorithm converged to a point of local infeasibility. Problem may be infeasible.
        USER_REQUESTED_STOP: The user call-back function intermediate_callback (see Section 3.3.4) returned false, i.e., the user code requested a premature termination of the optimization.
        DIVERGING_ITERATES: It seems that the iterates diverge.
        RESTORATION_FAILURE: Restoration phase failed, algorithm doesn't know how to proceed.
        ERROR_IN_STEP_COMPUTATION: An unrecoverable error occurred while IPOPT tried to compute the search direction.
        INVALID_NUMBER_DETECTED: Algorithm received an invalid number (such as NaN or Inf) from the NLP; see also option check_derivatives_for_naninf.
        INTERNAL_ERROR: An unknown internal error occurred. Please contact the IPOPT authors through the mailing list.
    n: (in), the number of variables in the problem (dimension of $ x$).
    x: (in), the final values for the primal variables, $ x_*$.
    z_L: (in), the final values for the lower bound multipliers, $ z^L_*$.
    z_U: (in), the final values for the upper bound multipliers, $ z^U_*$.
    m: (in), the number of constraints in the problem (dimension of $ g(x)$).
    g: (in), the final value of the constraint function values, $ g(x_*)$.
    lambda: (in), the final values of the constraint multipliers, $ \lambda_*$.
    obj_value: (in), the final value of the objective, $ f(x_*)$.
    ip_data and ip_cq are provided for expert users.
*/

virtual void finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i<n; i++) {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
}
  
  private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  HS071_NLP();
  HS071_NLP(const HS071_NLP&);
  HS071_NLP& operator=(const HS071_NLP&);
  //@}
  
};

#endif