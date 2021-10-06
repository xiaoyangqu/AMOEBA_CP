!---------------------------!
!   Code by: Qiantao Wang   !
!            Jul 2013       !
!---------------------------!
!
!   ==============================
!  |    MODULE chgpen             |
!  |                              |
!  |    This module calculates    |
!  |    the charge penetration    |
!  |    energy and force          |
!  |    between charge-charge     |
!  |    and charge-dipole         |
!  |    in AMOEBA                 |
!  |                              |
!  |    Two formulisms are        |
!  |    implemented:              |
!  |                              |
!  |    1) JP                     |
!  |    2) Truhlar's formular     |
!   ============================== 
! 
!
!  Module structure:
!
!     Public derived type
!        1) ecp
!
!

      module chgpen
        implicit none
     
        private
     
        ! subroutines and functions
        public :: ecp

        type ecp_type

           ! initialized in 'initprm'
           ! maxtyp long
           real*8  :: alp_jp(5000)  ! alpha in JP
           real*8  :: bet_jp(5000)  ! alpha in JP

           integer :: n_val(5000)   ! no. of valence electrons

           ! the flag controls the form of charge penetration
           ! 1 - JP
           ! 2 - Truhlar
           integer :: chgpen_form

           ! Default: true
           logical :: use_chgpen

           ! chgpen switching function
           real*8 :: switch_lo  ! lower bound
           real*8 :: switch_hi  ! higher bound

        end type ecp_type
     
        ! Derived type
        ! Public
        type(ecp_type) :: ecp

        ! Classes used in potential calculation
        public :: potcp

        type pot_type
           ! used in potential.x
           logical :: dochgpen

        end type pot_type

        ! Derived type
        ! Public
        type(pot_type) :: potcp

      end module chgpen
