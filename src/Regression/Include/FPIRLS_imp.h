#ifndef __FPIRLS_IMP_H__
#define __FPIRLS_IMP_H__

#include "FPIRLS.h"


/*********** FPIRLS_base Methods ************/

// Constructor
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim> & mesh, InputHandler & inputData, OptimizationData & optimizationData,  VectorXr mu0, bool scale_parameter_flag, Real scale_param):
  mesh_(mesh), inputData_(inputData), optimizationData_(optimizationData), regression_(inputData, optimizationData, mesh.num_nodes()), mu_(mu0), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param)
{
  //initialization of mu, current_J_values and past_J_values.  for(UInt j=0; j< inputData.getLambdaS().size() ; j++){
    //mu_.push_back(mu0);
    current_J_values = std::array<Real,2>{1,1};
    past_J_values = std::array<Real,2>{1,1};

};
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim> & mesh, const std::vector<Real> mesh_time, InputHandler & inputData, OptimizationData & optimizationData,  VectorXr mu0, bool scale_parameter_flag, Real scale_param):
  mesh_(mesh), mesh_time_(mesh_time), inputData_(inputData), optimizationData_(optimizationData), regression_(mesh_time, inputData, optimizationData, mesh.num_nodes()), mu_(mu0), scale_parameter_flag_(scale_parameter_flag), _scale_param(scale_param)
{
  //initialization of mu, current_J_values and past_J_values.  for(UInt j=0; j< inputData.getLambdaS().size() ; j++){
    //mu_.push_back(mu0);
    current_J_values = std::array<Real,2>{1,1};
    past_J_values = std::array<Real,2>{1,1};

};


// FPIRLS_base methods
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  //Initialize the containers size, as LambdaS_len
  const UInt LambdaS_len = (*optimizationData_.get_LambdaS_vector()).size();
  const UInt LambdaT_len = (*optimizationData_.get_LambdaT_vector()).size();
  
  // initialize the algorithm variables
  /*
  G_.resize(LambdaS_len);
  WeightsMatrix_.resize(LambdaS_len);
  pseudoObservations_.resize(LambdaS_len);
  */
  n_iterations.resize(LambdaS_len,LambdaT_len);
  _variance_estimates.resize(LambdaS_len,LambdaT_len);

  // Initialize the outputs. The temporal dimension is not implemented, for this reason the 2nd dimension is set to 1.
  //
  /* Quando hai fatto tutto qua devi adattare la seconda dimensione per il caso separabile, quindi con LambdaT:  <07-09-20, Elia Cunial> */
  if( this->inputData_.getCovariates()->rows() > 0 ) 
      _beta_hat.resize(LambdaS_len,LambdaT_len);
  _fn_hat.resize(LambdaS_len,LambdaT_len);
  _dof.resize(LambdaS_len,LambdaT_len);
  _solution.resize(LambdaS_len,LambdaT_len);
   
  _GCV.resize(LambdaS_len, LambdaT_len);//If GCV is not computed the vector stores -1
  _GCV.setConstant(-1);
  _J_minima.resize(LambdaS_len, LambdaT_len);
 

  if(isSpaceVarying)
  {
    FiniteElement< ORDER, mydim, ndim> fe;
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  for(UInt i=0 ; i < LambdaS_len ; i++){//for-cycle for each spatial penalization (lambdaS). 
      for(UInt j = 0; j < LambdaT_len; j++){
        current_J_values[0] = past_J_values[0] + 2*inputData_.get_treshold();
        current_J_values[1] = past_J_values[1] + 2*inputData_.get_treshold();
        this->optimizationData_.setCurrentLambda(i, j); // set right lambda for the current iteration.

  //  Rprintf("Start FPIRLS for the lambda number %d \n", i+1);

        n_iterations(i, j) = 0;
        while(stopping_criterion(i, j)){ 
        
          // STEP (1)
          
          compute_G();
          compute_Weights();
          compute_pseudoObs();

          // STEP (2)
          
          this->inputData_.updatePseudodata(pseudoObservations_, WeightsMatrix_);
          update_solution(i, j);

          // STEP (3)
          compute_mu(i, j);

          // Variance Estimate
          compute_variance_est(i, j);

          // update J
          past_J_values = current_J_values;
          current_J_values = compute_J(i, j);
          n_iterations(i, j)++;

        } //end while
        /*
        #ifdef R_VERSION_
        Rprintf("\t n. iterations: %d\n \n", n_iterations(i, j));
        #else
        std::cout<< "\t n. iterations: "<<n_iterations(i, j)<<"\n"<<std::endl;
        #endif
        _J_minima(i, j) = current_J_values[0]+current_J_values[1]; // compute the minimum value of the J fuctional 

        */
        if(this->optimizationData_.get_loss_function()=="GCV"){ // compute GCV if it is required
          compute_GCV(i, j);
        }

      }// end for

      }


}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::update_solution(UInt& lambdaS_index, UInt& lambdaT_index){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem. 
  regression_.recomputeWTW(); // at each iteration of FPIRLS W is updated, so WTW has to be recomputed as well.
  regression_.preapply(this->mesh_);
  regression_.apply();
  const SpMat * Psi = regression_.getpsi_(); // get Psi matrix. It is used for the computation of fn_hat.

  // get the solutions from the regression object.
  _solution(lambdaS_index, lambdaT_index) = regression_.getSolution()(0,0);
  _dof(lambdaS_index, lambdaT_index) = regression_.getDOF()(0,0);

  if(inputData_.getCovariates()->rows()>0){ 
    _beta_hat(lambdaS_index, lambdaT_index) = regression_.getBeta()(0,0);
  }

  _fn_hat(lambdaS_index, lambdaT_index) = (*Psi) *_solution(lambdaS_index,lambdaT_index).topRows(Psi->cols());

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_pseudoObs(){
  // compute pseudodata observations
  
  VectorXr first_addendum; // G_ii( z_i - mu_i)
  VectorXr g_mu; // g( mu_i )
  
  const VectorXr * z = inputData_.getInitialObservations();

  first_addendum.resize(mu_.size());
  g_mu.resize(mu_.size());

  //compute the vecotr first_addendum and g_mu
  for(auto i=0; i < mu_.size(); i++){
    g_mu(i) = link(mu_(i));
    first_addendum(i) = G_(i)*((*z)(i)-mu_(i));
  }

  pseudoObservations_ = first_addendum + g_mu;

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_G(){
  // compute the G matrix as G_ii = diag( g'(mu_i))
  
  G_.resize(mu_.size());

  for(UInt i = 0; i<mu_.size(); i++){
    G_(i) = link_deriv(mu_(i));
  }


}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_Weights(){
  // computed W elementwise (it is a diagonal matrix)

  WeightsMatrix_.resize( mu_.size());

  /* V(Y) = b''(theta) = 1/g'(mu) ==> W = 1/(g'(mu)^2*V(mu)) = 1/g'(mu):  <01-10-20, Elia Cunial> */
  for(auto i=0; i < mu_.size(); i++){
        WeightsMatrix_(i) = 1/G_(i);
        
  }

}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_mu(UInt& lambdaS_index, UInt& lambdaT_index){
  //compute mu as mu_i = g-1( w_ii*beta + fn_hat)

  VectorXr W_beta = VectorXr::Zero(mu_.size()); // initialize the vector w_ii*beta

  if(inputData_.getCovariates()->rows()>0)
    W_beta = (*(inputData_.getCovariates()))*_beta_hat(lambdaS_index, lambdaT_index);

  
  for(UInt j=0; j < W_beta.size(); j++){
    mu_(j) = inv_link(W_beta[j] + _fn_hat(lambdaS_index,lambdaT_index)(j));
  }
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::stopping_criterion(UInt& lambdaS_index, UInt& lambdaT_index){
  // return true if the f-PIRLS has to perform another iteration, false if it has to be stopped
  
  bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
  bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?
  bool do_stop_by_exception = false;

  if(n_iterations(lambdaS_index, lambdaT_index) > 0 && regression_.getDecInfo() != 0){
      do_stop_by_exception = true;
      std::cout << "For lambda_index=( " << lambdaS_index << ", " << lambdaT_index << ") eigen decomposition error: " << regression_.getDecInfo() << std::endl;
  }

  if(n_iterations(lambdaS_index, lambdaT_index) > inputData_.get_maxiter()){
    do_stop_by_iteration = true;
    std::cout << "For lambda_index=( " << lambdaS_index << ", " << lambdaT_index << ") max iterations reached " << std::endl;
  }

  if(n_iterations(lambdaS_index, lambdaT_index) > 1){
    if(abs(past_J_values[0]+past_J_values[1] - current_J_values[0] - current_J_values[1]) < inputData_.get_treshold()){
        do_stop_by_treshold = true;
   }
  }

  return !(do_stop_by_iteration || do_stop_by_treshold || do_stop_by_exception);
}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_J(UInt& lambdaS_index, UInt& lambdaT_index){
  // compute the functional J: it is divided in parametric and non parametric part
  Real parametric_value = 0;
  Real non_parametric_value = 0;
  Real tmp;

  VectorXr Lf;

  const VectorXr * z = inputData_.getInitialObservations();

  for(UInt i=0; i < mu_.size(); i++){
    tmp = sqrt( var_function( mu_(i)) ) * ((*z)(i) - mu_(i)) ;
    parametric_value += tmp*tmp;
  }

      Lf.resize(_solution(lambdaS_index,lambdaT_index).size()/2);
      for(UInt i=0; i< Lf.size(); i++){
        Lf(i) = _solution(lambdaS_index,lambdaT_index)(Lf.size() + i);
      }


  if(isSpaceVarying)
  {
      Lf = Lf - forcingTerm;
  }

  SpMat Int;
  std::array<Real, 2> lambdaST = { (*optimizationData_.get_LambdaS_vector())[lambdaS_index], (*optimizationData_.get_LambdaT_vector())[lambdaT_index] };
  if(inputData_.isSpaceTime()) { // && inputData_.getFlagParabolic()){

        UInt correction_size = inputData_.getFlagParabolic() ? - 1 : + 2;
        VectorXr intcoef(mesh_time_.size() + correction_size);
        intcoef.setConstant(mesh_time_[1]-mesh_time_[0]);
        intcoef(0) *= 0.5;
        SpMat IN(mesh_.num_nodes(), mesh_.num_nodes());
        IN.setIdentity();
        SpMat tmp = intcoef.asDiagonal().toDenseMatrix().sparseView();
        tmp = kroneckerProduct(tmp, IN);
        Int.resize(tmp.rows(), tmp.cols());
        Int = lambdaST[0]*(*regression_.getR0_())*tmp;
  }
  else{
      Int.resize(mesh_.num_nodes(), mesh_.num_nodes());
      Int = lambdaST[0]*(*regression_.getR0_());
  }
  non_parametric_value = Lf.transpose() * Int * Lf;

  std::array<Real,2> returnObject{parametric_value, non_parametric_value};

  return returnObject;
}


template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_GCV(UInt& lambdaS_index, UInt&lambdaT_index){

    if (optimizationData_.get_DOF_evaluation() != "not_required") //in this case surely we have already the dofs
    { // is DOF_matrix to be computed?
        regression_.computeDegreesOfFreedom(0, 0, (*optimizationData_.get_LambdaS_vector())[lambdaS_index], (*optimizationData_.get_LambdaT_vector())[lambdaT_index]);
        _dof(lambdaS_index,lambdaT_index) = regression_.getDOF()(0,0);
    }
    else 
        _dof(lambdaS_index,lambdaT_index) = regression_.getDOF()(lambdaS_index,lambdaT_index);

    const VectorXr * y = inputData_.getInitialObservations();
    Real GCV_value = 0;

    for(UInt j=0; j < y->size();j++)
        GCV_value += dev_function(mu_[j], (*y)[j]); //norm computation

    GCV_value *= y->size();

    GCV_value /= (y->size()-optimizationData_.get_tuning()*_dof(lambdaS_index,lambdaT_index))*(y->size()-optimizationData_.get_tuning()*_dof(lambdaS_index,lambdaT_index));

    _GCV(lambdaS_index, lambdaT_index) = GCV_value;

    //best lambda
    if(GCV_value < optimizationData_.get_best_value())
    {
    optimizationData_.set_best_lambda_S(lambdaS_index);
    optimizationData_.set_best_lambda_T(lambdaT_index);
    optimizationData_.set_best_value(GCV_value);
    }

}

template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::compute_variance_est(UInt & lambdaS_index, UInt & lambdaT_index){
  Real phi;
  if(this->scale_parameter_flag_ && this->optimizationData_.get_loss_function()!="GCV"){// if scale param should be
    _variance_estimates.resize(this->mu_.size(),0);
    const UInt n_obs = this->inputData_.getObservations()->size();

    //scale parameter computed as: mean((var.link(mu)*phi)/mu), and phi is computed as in Wood IGAM
      phi = (this->scale_parameter_flag_ )? this->current_J_values[0]/(n_obs - this->_dof(lambdaS_index,lambdaT_index) ) : _scale_param;
      for(UInt j=0; j < this->mu_.size(); j++){
        _variance_estimates( lambdaS_index, lambdaT_index ) += phi* this->var_function(this->mu_[j])/this->mu_[j];
      }
      _variance_estimates( lambdaS_index, lambdaT_index ) /= this->mu_.size();
  }else{
    _variance_estimates(lambdaS_index, lambdaT_index) = -1;
  }
}

/*********** FPIRLS apply template specialization ************/

// Laplace or Elliptic case
template <typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<InputHandler,ORDER, mydim, ndim>::apply(){

  FPIRLS_Base<InputHandler,ORDER, mydim, ndim>::apply(ForcingTerm());

}

// SpaceVarying case
template <UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<GAMDataEllipticSpaceVarying,ORDER, mydim, ndim>::apply(){

  this->isSpaceVarying = true;
  FPIRLS_Base<GAMDataEllipticSpaceVarying,ORDER, mydim, ndim>::apply(this->inputData_.getU());
}
#endif
