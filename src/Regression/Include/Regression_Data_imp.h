#ifndef __REGRESSION_DATA_IMP_H__
#define __REGRESSION_DATA_IMP_H__

// Laplace
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(
    SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
    SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
    SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold)
    : RegressionData(Rlocations, RbaryLocations, Robservations, Rorder,
                     Rcovariates, RBCIndices, RBCValues, RincidenceMatrix,
                     RarealDataAvg, Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
    threshold_ = REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    this->isGAM = true;
}

// PDE
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(
    SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
    SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices,
    SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
    SEXP Rmax_num_iteration, SEXP Rthreshold)
    : RegressionDataElliptic(Rlocations, RbaryLocations, Robservations, Rorder,
                             RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues,
                             RincidenceMatrix, RarealDataAvg, Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
    threshold_ = REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    this->isGAM = true;
}

// PDE SpaceVarying
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(
    SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
    SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices,
    SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
    SEXP Rmax_num_iteration, SEXP Rthreshold)
    : RegressionDataEllipticSpaceVarying(
          Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru,
          Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg,
          Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
    threshold_ = REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    this->isGAM = true;
}
//Laplace time
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(
    SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations,
    SEXP Robservations, SEXP Rorder, SEXP Rcovariates, SEXP RBCIndices,
    SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass,
    SEXP Rflag_parabolic, SEXP Ric, SEXP Rsearch, SEXP Rmax_num_iteration,
    SEXP Rthreshold, SEXP nochoice)
    : RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations,
                     Rorder, Rcovariates, RBCIndices, RBCValues,
                     RincidenceMatrix, RarealDataAvg, Rflag_mass,
                     Rflag_parabolic, Ric, Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
    threshold_ = REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    this->isGAM = true;
}

// PDE time
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(
    SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations,
    SEXP Robservations, SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc,
    SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
    SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,
    SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold)
    : RegressionDataElliptic(
          Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder,
          RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix,
          RarealDataAvg, Rflag_mass, Rflag_parabolic, Ric, Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
    threshold_ = REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    this->isGAM = true;
}

// PDE SpaceVarying time
template <typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(
    SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations,
    SEXP Robservations, SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
    SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix,
    SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,
    SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold)
    : RegressionDataEllipticSpaceVarying(
          Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder,
          RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues,
          RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Ric,
          Rsearch) {
    max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
    threshold_ = REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    this->isGAM = true;
}

#endif
