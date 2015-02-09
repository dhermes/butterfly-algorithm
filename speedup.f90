module speedup
  implicit none
  private
  public coeff_new_level, intermediate_coeffs, dft_kernel

contains

  subroutine dft_kernel(t, s, K)
    ! NOTE: This assumes dp == 8 --> complex is two times that.
    real*8, intent(in) :: t, s
    complex*16, intent(out) :: K

    K = exp(cmplx(0.0d0, -1.0d0, 16) * t * s)

  end subroutine dft_kernel

  subroutine intermediate_coeffs(D_tau_sigma, D_tau_sigma_prime, tau, &
                                 tau_plus, sigma, sigma_prime, alpha, &
                                 factorial_values, R, &
                                 sub_result, sub_result_prime)
    ! Compute D(tau_+, sigma, alpha) and D(tau_+, sigma', alpha)
    ! 2-space.

    integer, intent(in) :: alpha, R
    complex*16, intent(in), dimension(0:R - 1) :: D_tau_sigma
    complex*16, intent(in), dimension(0:R - 1) :: D_tau_sigma_prime
    real*8, intent(in) :: tau, tau_plus, sigma, sigma_prime
    real*8, intent(in), dimension(0:R - 1) :: factorial_values
    complex*16, intent(out) :: sub_result, sub_result_prime

    ! Variables outside of signature.
    integer :: gamma
    real*8 :: comb_val
    complex*16 :: K

    sub_result = cmplx(0.0d0, 0.0d0, 16)
    sub_result_prime = cmplx(0.0d0, 0.0d0, 16)

    ! gamma == beta + alpha ==> alpha <= gamma < R
    do gamma = alpha, R - 1
        comb_val = ( &
            (factorial_values(gamma) * (tau_plus - tau)**(gamma - alpha)) / &
            (factorial_values(gamma - alpha) * factorial_values(alpha)))
        sub_result = sub_result + comb_val * D_tau_sigma(gamma)
        sub_result_prime = (sub_result_prime + &
                            comb_val * D_tau_sigma_prime(gamma))
    enddo

    call dft_kernel(tau_plus - tau, sigma, K)
    sub_result = sub_result * K
    call dft_kernel(tau_plus - tau, sigma_prime, K)
    sub_result_prime = sub_result_prime * K

  end subroutine intermediate_coeffs

  subroutine coeff_new_level(D_tau_sigma, D_tau_sigma_prime, tau, &
                             tau_plus, sigma, sigma_prime, sigma_minus, &
                             alpha, factorial_values, R, new_coeff)
    integer, intent(in) :: alpha, R
    complex*16, intent(in), dimension(0:R - 1) :: D_tau_sigma
    complex*16, intent(in), dimension(0:R - 1) :: D_tau_sigma_prime
    real*8, intent(in) :: tau, tau_plus, sigma, sigma_prime, sigma_minus
    real*8, intent(in), dimension(0:R - 1) :: factorial_values
    complex*16, intent(out) :: new_coeff

    new_coeff = cmplx(0.0d0, 0.0d0, 16)

  end subroutine coeff_new_level

end module speedup
