
!This is a test for your code.
!
!Just type these commands,
!
! ifort test5.f90 -mkl -$delm 
! ./a.out

module useful_tools

use de_model_init
use de_chisqs
implicit none



contains

  !------------------------------------------
  ! get the value of rho_de(z)/rho_de(z=0) 
  !  for maurice hde
  !------------------------------------------   
       DOUBLE PRECISION FUNCTION de_get_rhoz(inputa)
                DOUBLE PRECISION, intent(in) :: inputa
                DOUBLE PRECISION :: ez, z, a, ezsq, desq, othersq
                if(inputa .le. de_minintpla) then
                        a = de_minintpla
                else
                        a = inputa
                endif
                z = 1.0/a - 1.d0
                ez = de_get_ez(z)
                othersq = de_CP%Om0 * a**(-3.0d0) + de_CP%Or0 * a**(-4.0d0)+ de_CP%Ok0 * a**(-2.0d0)
                desq = ez**2.0 - othersq
                de_get_rhoz = desq/(1.0-de_CP%Om0-de_CP%Or0-de_CP%Ok0) !de_CP%Ode0
!                       de_get_rhoz = de_get_rhoz * (de_minintpla/inputa)**3.0
       END FUNCTION de_get_rhoz

end module useful_tools

program main

use de_model_init
use de_chisqs
use useful_tools

implicit none

	double precision :: y, y2
 	double precision :: omegam, z, a, ez, qz, q_ez_fun, H_residual
 	integer :: i, iz
 	character(len=1000) :: tmpstr
 
 	pr_chisq_info = .false.
 

!==========================================
! de_mauricehde
!==========================================

	print *
	print *
	
!        DO i = 1, 1
	de_model_lab = de_mauricehde_lab
	de_CP%Ob0hsq    =  0.02253
	omegam 		=  0.30 !0.26+i*0.002!0.284936E+00
	!de_CP%h		=  0.711833E+00
	de_CP%h		=  0.70
	de_CP%alpha	=  0.142125E+01
	de_CP%beta	=  0.325121E+01  
	de_CP%Odm0 	=  omegam - de_CP%Ob0hsq/de_CP%h**2.0
	de_CP%Ok0	= 0.0

	call de_init()
	write(tmpstr,'(A,f4.2,A,f4.2,A,f10.8)') '### omegam = ', real(omegam), '; h = ', de_CP%h, &
	  '; omega_radiaton = ', de_CP%Or0

!	y = de_chisq_g06(.false., .true.)
!	y = de_chisq_jla()
!	y = de_chisq_snls3() !+ de_chisq_sdssdr7_old() + de_chisq_wmap7() + de_chisq_h_Riess()
!	y = de_chisq_union2p1() !+ de_chisq_wmap7()
!	y = de_chisq_wmap7()
!	y = de_chisq_planck()

	open(unit=987,file='z_ez_qz__mauricehde.txt')
	write(987,*) trim(adjustl(tmpstr))
	do iz = 1, de_num_intpl, 20
		z=de_zi(iz)
		!z = ia * 0.1d0
		a=1.0/(1.0+z)
		ez=1.0/de_inv_e(z)
		qz=de_mauricehde_q(a,ez)
!		q_ez_fun=(1-qz)*ez*ez
!		H_residual=de_CP%Om0/a**3 + de_CP%Or0/a**4 + (1-qz)*ez*ez / 3.0 - ez*ez
		!print *, 'z / h / q = ', real(z), real(ez), real(de_mauricehde_q(a,ez))
		write(987,*)  real(z), real(ez), real(qz)!, real(q_ez_fun), real(H_residual)
	enddo
	close(987)
	
	! Compute quantities for Lambda CDM
	de_model_lab = de_lcdm_lab
	call de_init()
!	y2 = de_chisq_snls3() !+ de_chisq_sdssdr7_old() + de_chisq_wmap7() + de_chisq_h_Riess()
!	y2 = de_chisq_union2p1() !+ de_chisq_wmap7()
!	y2 = de_chisq_wmap7()
!	y2 = de_chisq_planck()

	open(unit=987,file='z_ez_qz__lcdm.txt')
	write(987,*) trim(adjustl(tmpstr))
	do iz = 1, de_num_intpl, 20
		z=de_zi(iz)
		!z = ia * 0.1d0
		a=1.0/(1.0+z)
		ez=1.0/de_inv_e(z)
		qz= -1.0 - (a/ez) * (-3*de_CP%Om0*a**(-4) -4*de_CP%Or0*a**(-5))/(2*ez)
!		q_ez_fun=(1-qz)*ez*ez
		!print *, 'z / h / q = ', real(z), real(ez), real(qz)
		write(987,*)  real(z),  real(ez), real(qz)!,  real(q_ez_fun)
	enddo
	close(987)

!	print *, 'omegam /chisqs of HDE/LambdaCDM =', real(omegam), real(y), real(y2)
!	ENDDO
!	write(*,*) ""
!	write(*,*) "====================================="
!	write(*,*) "  Resut of Maurice's HDE:"
!	write(*,*) "     Total chisq (lnlike) = ", y, y/2.0
!	write(*,*) "     Expected chisq = ", 424.911450626633d0
!	write(*,*) "====================================="	

end program main
