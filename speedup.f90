module speedup
  use types, only: dp
  implicit none
  private
  public coeff_new_level, intermediate_coeffs

contains

  subroutine intermediate_coeffs(D_tau_sigma, D_tau_sigma_prime, tau, &
                                 tau_plus, sigma, sigma_prime, alpha, R, &
                                 sub_result, sub_result_prime)
    ! Compute D(tau_+, sigma, alpha) and D(tau_+, sigma', alpha)
    ! 2-space.

    integer, intent(in) :: alpha, R
    complex(dp), intent(in) :: D_tau_sigma(R)
    complex(dp), intent(in) :: D_tau_sigma_prime(R)
    real(dp), intent(in) :: tau, tau_plus, sigma, sigma_prime
    complex(dp), intent(out) :: sub_result, sub_result_prime

    sub_result = 0.0_dp
    sub_result_prime = 0.0_dp

  end subroutine intermediate_coeffs

  subroutine coeff_new_level(D_tau_sigma, D_tau_sigma_prime, tau, &
                             tau_plus, sigma, sigma_prime, sigma_minus, &
                             alpha, R, new_coeff)
    integer, intent(in) :: alpha, R
    complex(dp), intent(in) :: D_tau_sigma(R)
    complex(dp), intent(in) :: D_tau_sigma_prime(R)
    real(dp), intent(in) :: tau, tau_plus, sigma, sigma_prime, sigma_minus
    complex(dp), intent(out) :: new_coeff

    new_coeff = 0.0_dp

  end subroutine coeff_new_level

end module speedup
