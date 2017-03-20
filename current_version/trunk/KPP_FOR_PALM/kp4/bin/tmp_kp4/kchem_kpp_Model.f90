MODULE kchem_kpp_Model

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Completely defines the model kchem_kpp
!    by using all the associated modules
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE kchem_kpp_Precision
  USE kchem_kpp_Parameters
  USE kchem_kpp_Global
  USE kchem_kpp_Function
  USE kchem_kpp_Integrator
  USE kchem_kpp_Rates
  USE kchem_kpp_Jacobian
  USE kchem_kpp_Hessian
  USE kchem_kpp_Stoichiom
  USE kchem_kpp_LinearAlgebra
  USE kchem_kpp_Monitor
  USE kchem_kpp_Util

END MODULE kchem_kpp_Model

