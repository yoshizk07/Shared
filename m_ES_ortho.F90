!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 603 $)
!
!  MODULE: m_Electronic_Structure
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004, April/15/2006, September/02/2008
!  FURTHER MODIFICATION: T. Yamasaki, T. Uda and T. Ohno, September 2009 (MGS_DGEMM)
!  FURTHER MODIFICATION: T. Yamasaki and T. Yamamoto,   October 2009  (NONLOCAL_DGEMM)
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!

!   This module has been revised for the GAMMA point (k=(0,0,0)) by T. Yamasaki
!  in April 2006. Number of operations for the Gamma point have been tremendously
!  reduced in subroutines of m_ES_betar_dot_Wfs, m_ES_Vnonlocal_W, and
!  m_ES_modified_gram_schmidt.
!
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_DO__
#   define __TIMER_DO_START(a)   call timer_sta(a)
#   define __TIMER_DO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_DO_START(a)
#   define __TIMER_DO_STOP(a)
#endif
#ifdef __TIMER_COMM__
#   define __TIMER_COMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif
#ifdef __TIMER_DGEMM__
#   define __TIMER_DGEMM_START(a)  call timer_sta(a)
#   define __TIMER_DGEMM_STOP(a)   call timer_end(a)
#else
#   define __TIMER_DGEMM_START(a)
#   define __TIMER_DGEMM_STOP(a)
#endif

!  2010.03
!  the performance of the Gram-Schmidt orthonormalization procedure
!  has been tuned through a joint effort with RIKEN.

#ifndef SX
#define DGEMM__       DGEMM
#endif

#ifndef NO_MGS_DGEMM
#define MGS_DGEMM
#endif

module m_ES_ortho
  use m_Control_Parameters, only : nspin,ipri,kimg, neg, printable, af &
       &                         , m_CtrlP_cachesize &
       &                         , sw_gep &
#ifdef SAVE_FFT_TIMES
       &                         , nblocksize_mgs, nblocksize_mgs_is_given, sw_save_fft
#else
       &                         , nblocksize_mgs, nblocksize_mgs_is_given
#endif
  use m_Const_Parameters,   only : DP, EXECUT, ON, OFF, DELTA, PAI2 &
       &                         , ORTHOGONALIZATION, ORTHONORMALIZATION, NORMALIZATION &
       &                         , NORMCONSERVATION, VDB, VANDERBILT_TYPE &
       &                         , OTHER_BANDS, SAME_BAND, GAMMA, OLD
  use m_Electronic_Structure,only: zaj_l, neordr, nrvf_ordr, fsr_l, fsi_l &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
       &                         , nblocksize_mgs_default
  use m_Files,              only : nfout
  use m_Kpoints,            only : kv3,k_symmetry
  use m_Parallelization,    only : MPI_CommGroup,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,ista_k,iend_k &
       &                         , ierr,mp_e,nis_e,nie_e,nel_e &
       &                         , ista_g1k,iend_g1k,np_g1k,mp_g1k,nis_g1k,nel_g1k &
       &                         , np_fs, mp_fs, ista_fs, iend_fs, nis_fs, nie_fs, nel_fs  &
       &                         , myrank_g, nrank_g &
       &                         , is_ffth, ista_ffth, nel_ffth, np_ffth, mp_ffth, neg_g_all
  use m_PlaneWaveBasisSet,  only : kg1,kg,iba
  use m_PseudoPotential,    only : nlmt,nlmta,lmta,lmtt,ltp,q &
       &                         , modnrm,nac,fqwei,nlmta1,nlmta2 &
       &                         , nac_p, fqwei_p, nlmta1_p, nlmta2_p
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Electronic_Structure, only : fsr_l_2D, fsi_l_2D           &
 &                                 , fsr_ball, fsi_ball, zaj_ball           &
 &                                 , m_ES_wd_zaj_small_portion_3D
  use m_ES_nonlocal,          only : m_ES_betar_dot_WFs_4_each_k_3D
  use m_Parallelization,      only : lrank, nbsn, nbsn_sta, nbsn_end, neg_g           &
 &                                 , nbs_num, nbsn_num, nbs_sta, nbs_end              &
 &                                 , mpi_kg_world, mpi_ke_world           &
 &                                 , map_k                   &
 &                                 , mp_kngp, is_kngp, ie_kngp, np_kngp, nel_kngp     &
 &                                 , ista_spin, iend_spin

! ================================ added by K. Tagami ================ 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor
  use m_PseudoPotential,      only : fqwei_noncl, fqwei_p_noncl
! ==================================================================== 11.0
  use mpi

  implicit none

  integer, private, parameter                         :: sw_timing_2ndlevel = ON

  integer                                             :: np_g1k_x     ! np_g1k_x = np_g1k(ik) | np_g1k(ik)+1
  real(kind=DP),private,allocatable,target,dimension(:,:,:)  :: psi_t  ! d(np_g1k_x,neg,kimg) a work_array for mgs
  real(kind=DP),private,allocatable,dimension(:,:)    :: p1Sp2
  integer                                             :: np_fs_x      ! np_fs_x = np_fs | np_fs+1
  real(kind=DP),private,allocatable,dimension(:,:)    :: bpr_t, bpi_t ! d(np_fs_x,neg) work arrays for mgs
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
  real(kind=DP),private,allocatable,dimension(:)      :: psi_ii,psi_ir
  real(kind=DP),private,allocatable,dimension(:)      :: bp_ii,bp_ir
#endif
#ifdef MGS_DGEMM
!-$ tune FUJITSU for block ----->>
  real(kind=DP),private,allocatable,dimension(:,:)    :: bpr_tw1, bpr_tw2, bpi_tw1, bpi_tw2 ! d(nac_p,neg) work arrays for mgs
  integer,private                                     :: NB
#ifdef SX
  integer,private                                     :: MB, LMB
#endif
! <--
!-$ tune FUJITSU for block <<-----
#endif
! <--

!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

! ============================= added by K. Tagami ================ 11.0
  real(kind=DP), private, allocatable, target :: psi_t_noncl(:,:,:,:)
  real(kind=DP), private, allocatable :: bpr_t_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpi_t_noncl(:,:,:)

#ifdef MGS_DGEMM
  real(kind=DP), private, allocatable :: bpr_tw1_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpr_tw2_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpi_tw1_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpi_tw2_noncl(:,:,:)
#endif
! ================================================================= 11.0

contains
!  1-18. m_ESortho_mgs_alloc
!        - mgs_vdb_alloc, - mgs_nrc_alloc
!  1-19. m_ESortho_mgs_dealloc
!        - mgs_vdb_dealloc, - mgs_nrc_dealloc
!  2-28. m_ES_modified_gram_schmidt
!  2-29. m_ES_orthogonalize_SD_to_WFs -> (2-30)
!  2-30. mgs_sd2wf_each_k_G           <- (2-29)  -> (2-34), (2-44)
!        - broadcast_fs, - Psi1SPhi2_t, - modify_bsd_and_phi_t, - alloc_phi_w_and_brd_w
!        - dealloc_phi_w_and_brd_w, - WSW, - normalize_bsd_and_phi, - Psi1SPhi2
!        - modify_bsd_and_phi
!  2-31. orthogonalize_SD             -> (2-33)
!  2-32. m_ES_MGS_4_each_k            -> (2-17), (2-33)
!        - wd_title_of_the_operation
!  2-33. mgs_4_each_k_G               <- (2-31), (2-32)  -> (2-34), (2-44)
!        - WSW_t, - normalize_bp_and_psi_t, - W1SW2_t_r
!        - modify_bp_and_psi_t_r, - substitute_jto_ib2back, - W1SW2_t
!        - modify_bp_and_psi_t, - alloc_and_brd_w, - dealloc_and_brd_w
!        - WSW, - normalize_bp_and_psi, - W1SW2, - modify_bp_and_psi
!  2-34. set_npzri                    <- (2-30), (2-33)
!  2-37. m_ES_W_transpose_back2
!  2-38. m_ES_W_transpose_r
!  2-39. m_ES_W_transpose_back_r
!  2-42. m_ES_F_transpose_r
!  2-43. m_ES_F_transpose_back_r
!   2-44. m_ES_F_transpose_r_3D
!   2-45. m_ES_F_transpose_back_r_3D
!   2-46. m_ES_W_transpose_r_3D
!   2-47. m_ES_W_transpose_back_r_3D

!  2-44. cp_bp_and_psi_2_brd_w2      <- (2-30), (2-33)

#ifdef VPP
#define _ODD_BOUNDARY_
#endif
#ifdef SX
#define _ODD_BOUNDARY_
#endif
#ifdef NEC_TUNE1
#define _ODD_BOUNDARY_
#endif


  subroutine m_ESortho_set_np_g1k_x(ik)
    integer, intent(in) :: ik
!   OUTPUT : np_g1k_x
!
!!$    np_g1k_x = np_g1k(ik) + 1
    np_g1k_x = np_g1k(ik)
  end subroutine m_ESortho_set_np_g1k_x

  subroutine m_ESortho_set_np_fs_x()
!!$    np_fs_x = np_fs + 1
    np_fs_x = np_fs
  end subroutine m_ESortho_set_np_fs_x

  subroutine m_ESortho_mgs_alloc(ik)
    integer, intent(in) :: ik
    call m_ESortho_set_np_g1k_x(ik) ! -> np_g1k_x
!!$    call set_np_g1k_x() ! -> np_g1k_x

    if(modnrm == EXECUT) then
       call m_ESortho_set_np_fs_x() ! -> np_fs_x
!!$       call set_np_fs_x() ! -> np_fs_x
! ============================== modified by K. Tagami ============= 11.0
!       call mgs_vdb_alloc
       call mgs_vdb_alloc(np_e)
! =================================================================== 11.0
    else
! ============================== modified by K. Tagami ============= 11.0
!       call mgs_nrc_alloc
       call mgs_nrc_alloc(np_e)
! ===================================================================== 11.0
    end if
  contains
!!$    subroutine set_np_g1k_x()
!!$!   OUTPUT : np_g1k_x
!!$!
!!$!BRANCH_P ORG_Parallel
!!$#ifdef _ODD_BOUNDARY_
!!$      if(mod(np_g1k(ik),2) == 0) then
!!$         np_g1k_x = np_g1k(ik) + 1
!!$      else
!!$         np_g1k_x = np_g1k(ik)
!!$      end if
!!$#else
!!$!BRANCH_P_END ORG_Parallel
!!$      np_g1k_x = np_g1k(ik) + 1
!!$!BRANCH_P ORG_Parallel
!!$#endif
!!$!BRANCH_P_END ORG_Parallel
!!$    end subroutine set_np_g1k_x
!!$
!!$    subroutine set_np_fs_x()
!!$
!!$!BRANCH_P ORG_Parallel
!!$#ifdef _ODD_BOUNDARY_
!!$      if(mod(np_fs,2) == 0) then
!!$         np_fs_x = np_fs+1
!!$      else
!!$         np_fs_x = np_fs
!!$      end if
!!$#else
!!$!BRANCH_P_END ORG_Parallel
!!$      np_fs_x = np_fs + 1
!!$!BRANCH_P ORG_Parallel
!!$#endif
!!$!BRANCH_P_END ORG_Parallel
!!$    end subroutine set_np_fs_x

! ============================== added by K. Tagami ================== 11.0
    subroutine mgs_vdb_alloc(neg)
      integer,intent(in) :: neg
      integer :: kimg_t

      allocate(psi_t(np_g1k_x,neg,kimg))

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      allocate(psi_ir(np_g1k_x))
      if(kimg == 2) allocate(psi_ii(np_g1k_x))
      if(np_g1k(ik) < np_g1k_x) then
         psi_ir(np_g1k_x) = 0.d0
         if(kimg==2) psi_ii(np_g1k_x) = 0.d0
      end if
#endif

      if((k_symmetry(ik) == GAMMA .and. kimg == 2) .or. kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(p1Sp2(neg,kimg_t))

      if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
         kimg_t = 1
      else
         kimg_t = 2
      end if


      if(kimg_t == 1) then
         allocate(bpr_t(np_fs_x,neg))
         allocate(bpi_t(1,1))
#ifdef MGS_DGEMM
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         if(nac_p > 0) then
!  =============================================================================
         allocate(bpr_tw1(nac_p,neg))
         allocate(bpr_tw2(nac_p,neg))
         allocate(bpi_tw1(nac_p,neg))
         allocate(bpi_tw2(nac_p,neg))
!         allocate(bpi_tw1(1,1))
!         allocate(bpi_tw2(1,1))
         bpr_tw1 = 0.0
         bpr_tw2 = 0.0
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         end if
!  =============================================================================
#endif
      else
         allocate(bpr_t(np_fs_x,neg))
         allocate(bpi_t(np_fs_x,neg))
#ifdef MGS_DGEMM
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         if(nac_p > 0) then
!  =============================================================================
         allocate(bpr_tw1(nac_p,neg))
         allocate(bpr_tw2(nac_p,neg))
         allocate(bpi_tw1(nac_p,neg))
         allocate(bpi_tw2(nac_p,neg))
         bpr_tw1 = 0.0
         bpr_tw2 = 0.0
         bpi_tw1 = 0.0
         bpi_tw2 = 0.0
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         end if
!  =============================================================================
#endif
      end if

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      allocate(bp_ir(max(np_fs,np_g1k_x)));bp_ir=0.d0
      if(kimg_t == 2) then
         allocate(bp_ii(max(np_fs,np_g1k_x)));bp_ii=0.d0
      endif
!!$      allocate(bp_ii(np_g1k_x))
      if(np_g1k(ik) < np_g1k_x) then
         bp_ir(np_g1k(ik)+1:np_g1k_x) = 0.d0
!!$         if(kimg==2) bp_ii(np_g1k_x) = 0.d0
         if(kimg_t == 2) bp_ii(np_g1k(ik)+1:np_g1k_x) = 0.d0
!!$         bp_ii(np_g1k(ik)+1:np_g1k_x) = 0.d0
      end if
#endif
    end subroutine mgs_vdb_alloc


    subroutine mgs_nrc_alloc(neg)
      integer,intent(in) :: neg
      integer :: kimg_t
      allocate(psi_t(np_g1k_x,neg,kimg))
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      allocate(psi_ir(np_g1k_x))
      if(kimg==2) allocate(psi_ii(np_g1k_x))
      if(np_g1k(ik) < np_g1k_x) then
         psi_ir(np_g1k_x) = 0.d0
         if(kimg==2) psi_ii(np_g1k_x) = 0.d0
      end if
#endif

      if((k_symmetry(ik) == GAMMA .and. kimg == 2).or.kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(p1Sp2(neg,kimg_t))
    end subroutine mgs_nrc_alloc
  end subroutine m_ESortho_mgs_alloc


  subroutine m_ESortho_mgs_dealloc()
    if(modnrm == EXECUT) then
! ================================ modified by K. Tagami ========== 11.0
!!       call mgs_vdb_dealloc()
!
       if ( noncol ) then
         call mgs_vdb_dealloc_noncl()
       else
         call mgs_vdb_dealloc()
       endif
! ================================================================= 11.0
    else
! ================================ modified by K. Tagami ========== 11.0
!       call mgs_nrc_dealloc()

       if ( noncol ) then
          call mgs_nrc_dealloc_noncl()
       else
          call mgs_nrc_dealloc()
       endif
! ==================================================================== 11.0
    end if
  contains
    subroutine mgs_vdb_dealloc()
      deallocate(bpr_t)
#ifdef MGS_DGEMM
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
!     deallocate(bpr_tw1)
!     deallocate(bpr_tw2)
      if(allocated(bpr_tw1)) deallocate(bpr_tw1)
      if(allocated(bpr_tw2)) deallocate(bpr_tw2)
!  =============================================================================
#endif
      if(allocated(bpi_t)) deallocate(bpi_t)
#ifdef MGS_DGEMM
      if(allocated(bpi_tw1)) deallocate(bpi_tw1)
      if(allocated(bpi_tw2)) deallocate(bpi_tw2)
#endif
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      deallocate(psi_ir)
      if(kimg == 2) deallocate(psi_ii)
#endif
      deallocate(psi_t)
      deallocate(p1Sp2)
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      if(allocated(bp_ir)) deallocate(bp_ir)
      if(allocated(bp_ii)) deallocate(bp_ii)
#endif
    end subroutine mgs_vdb_dealloc

! ============================= added by K. Tagami ================= 11.0
    subroutine mgs_vdb_dealloc_noncl()
      deallocate(bpr_t_noncl)
#ifdef MGS_DGEMM
      deallocate(bpr_tw1_noncl)
      deallocate(bpr_tw2_noncl)
#endif
      if(allocated(bpi_t_noncl)) deallocate(bpi_t_noncl)
#ifdef MGS_DGEMM
      if(allocated(bpi_tw1_noncl)) deallocate(bpi_tw1_noncl)
      if(allocated(bpi_tw2_noncl)) deallocate(bpi_tw2_noncl)
#endif
      deallocate(psi_t_noncl)
      deallocate(p1Sp2)
!!$#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!!$      if(mod_pot == VANDERBILT_TYPE) then
!!$         if(allocated(bp_ir)) deallocate(bp_ir)
!!$         if(allocated(bp_ii)) deallocate(bp_ii)
!!$      end if
!!$#endif
    end subroutine mgs_vdb_dealloc_noncl
! ================================================================= 11.0

    subroutine mgs_nrc_dealloc()
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      deallocate(psi_ir)
      if(kimg == 2) deallocate(psi_ii)
#endif
      deallocate(psi_t)
      deallocate(p1Sp2)
    end subroutine mgs_nrc_dealloc

! ============================== added by K. Tagami =============== 11.0
    subroutine mgs_nrc_dealloc_noncl()
      deallocate(psi_t_noncl)
      deallocate(p1Sp2)
    end subroutine mgs_nrc_dealloc_noncl
! ================================================================== 11.0

  end subroutine m_ESortho_mgs_dealloc

  subroutine m_ES_modified_gram_schmidt(nfout)
  use m_IterationNumbers,     only : iteration
    integer, intent(in) :: nfout
    integer :: ik,is

    do is = ista_spin, iend_spin, (af+1)
!    do ik = ista_k, iend_k, af+1                            ! MPI
    do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle          ! MPI
       if(ipri>=2) call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- before GS --",16)

       if(sw_gep == OFF .or. iteration<2)then
          call m_ES_MGS_4_each_k(nfout,ik,mode=ORTHONORMALIZATION)
       else
          call m_ES_MGS_4_each_k(nfout,ik,mode=NORMALIZATION)
       endif
       !    ~~~~~~~~~~~~~~~~~~~~~~
       if(ipri>=2) call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- after GS --",15)

    end do
    end do

  end subroutine m_ES_modified_gram_schmidt

!!$!BRANCH_P 3D_Parallel
!!$  subroutine m_ES_modified_gram_schmidt_3D(nfout)
!!$    integer, intent(in) :: nfout
!!$    integer :: ik
!!$
!!$    do ik = ista_k, iend_k, af+1                            ! MPI
!!$       if(ipri>=2) call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- before GS --",16)
!!$
!!$       call m_ES_MGS_4_each_k(nfout,ik,mode=ORTHONORMALIZATION)
!!$       !    ~~~~~~~~~~~~~~~~~~~~~~
!!$       if(ipri>=2) call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- after GS --",15)
!!$    end do
!!$
!!$  end subroutine m_ES_modified_gram_schmidt_3D
!!$!BRANCH_P_END 3D_Parallel

  subroutine m_ES_orthogonal_phi_to_WFs(ik,wfsd_l,bsdr_l,bsdi_l)
    integer, intent(in)                                     :: ik
    real(kind=DP),intent(inout),dimension(kg1,np_e,ik:ik,kimg) :: wfsd_l
    real(kind=DP),intent(inout),dimension(np_e,nlmta,ik:ik)    :: bsdr_l,bsdi_l

    integer  :: id_sname = -1, id_sname2 = -1

    call tstatc0_begin('m_ES_orthogonal_phi_to_WFs ',id_sname,1)

    call tstatc0_begin('mpi_barrier(ortho_phi_to_WFs) ',id_sname2)
    call mpi_barrier(mpi_k_world(myrank_k),ierr)
    call tstatc0_end(id_sname2)

    if(modnrm == EXECUT) then
       call mgs_phi2wf_each_k_G(ik,wfsd_l,ORTHOGONALIZATION,bsdr_l,bsdi_l,mod_pot=VDB)
    else
       call mgs_phi2wf_each_k_G(ik,wfsd_l,ORTHOGONALIZATION,mod_pot=NORMCONSERVATION)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_ES_orthogonal_phi_to_WFs

  subroutine m_ES_orthogonalize_SD_to_WFs_3D(ik,to_which_band,wfsd_l,bsdr_l,bsdi_l)
    integer, intent(in)                                     :: ik,to_which_band
    real(kind=DP),intent(inout),dimension(maxval(np_g1k),np_e,ik:ik,kimg) :: wfsd_l
    real(kind=DP),intent(inout),dimension(np_e,nlmta,ik:ik)    :: bsdr_l,bsdi_l

    integer  :: id_sname = -1
    call tstatc0_begin('m_ES_orthogonalize_SD_to_WFs_3D ',id_sname)

!!$    call phi_alloc(ik)
!!$    call m_ES_mgs_alloc(ik)
    if(modnrm == EXECUT) then
       call mgs_sd2wf_each_k_G_3D(ik,to_which_band,wfsd_l,bsdr_l,bsdi_l,mod_pot=VDB)
    else
       call mgs_sd2wf_each_k_G_3D(ik,to_which_band,wfsd_l,mod_pot=NORMCONSERVATION)
    end if
    ! -(m_E.S.)   ->wfsd_l,bsdr_l,bsdi_l!!$    if(ik == iend_k) call phi_dealloc()    !-(m_E.S.)
!!$    call m_ES_mgs_dealloc

!!$    if(ik == iend_k) call phi_dealloc()    !-(m_E.S.)

    call tstatc0_end(id_sname)
  end subroutine m_ES_orthogonalize_SD_to_WFs_3D


  subroutine mgs_sd2wf_each_k_G_3D(ik,to,phi_l,bsdr_l,bsdi_l,mod_pot)
    integer, intent(in)             :: ik,to
    real(kind=DP), intent(inout)    :: phi_l(maxval(np_g1k),np_e,ik:ik,kimg)
    real(kind=DP), optional,intent(inout), dimension(np_e,nlmta,ik:ik) :: bsdr_l,bsdi_l
    integer, intent(in)             :: mod_pot
    integer               :: i

    integer :: kimg_t
    real(kind=DP),allocatable,dimension(:,:)    :: bp_w
#ifndef ZAJ_BALL_ALLREDUCE
! ==============================================================================
    real(kind=DP), allocatable, dimension(:,:,:)   :: wk_zaj !d(max_g1k,mp_e,kimg)
!!$    real(kind=DP), allocatable, dimension(:,:,:,:) :: wk_mpi

    integer :: max_g1k
    integer :: ierr, nb,kb, jb, ib

    if(nrank_e == 1) then
       zaj_ball(:,:,ik,:) = zaj_l(:,:,ik,:)
    else
       max_g1k = maxval(np_g1k(:))
       allocate(wk_zaj(max_g1k,mp_e,kimg),stat=ierr)
       do nb = 0, nrank_e-1
          if(nb == myrank_e) then
             do kb = 1, kimg
                do jb = 1, np_e
                   wk_zaj(:,jb,kb) = zaj_l(:,jb,ik,kb)
                end do
             end do
             if(np_e+1 < mp_e) wk_zaj(:,np_e+1:mp_e,:) = 0.d0
          end if
          call mpi_bcast(wk_zaj,max_g1k*mp_e*kimg,mpi_double_precision,nb,mpi_kg_world,ierr)
          do kb = 1, kimg
             do jb = nis_e(nb), nie_e(nb)
                do ib = 1, np_g1k(ik)
                   zaj_ball(ib,jb,ik,kb) = wk_zaj(ib,jb-nis_e(nb)+1,kb)
                end do
             end do
          end do
       end do
       deallocate(wk_zaj)
!!$       max_g1k = maxval(np_g1k(:))
!!$       allocate(wk_zaj(max_g1k,mp_e,kimg),stat=ierr); wk_zaj = 0.d0
!!$       allocate(wk_mpi(max_g1k,mp_e,kimg,0:nrank_e-1),stat=ierr)
!!$       do kb = 1, kimg
!!$          do jb = 1, np_e
!!$             do ib = 1, np_g1k(ik)
!!$                wk_zaj(ib,jb,kb) = zaj_l(ib,jb,ik,kb)
!!$             end do
!!$          end do
!!$       end do
!!$       call mpi_allgather(wk_zaj, max_g1k*mp_e*kimg, MPI_DOUBLE_PRECISION &
!!$            &            ,wk_mpi, max_g1k*mp_e*kimg, MPI_DOUBLE_PRECISION, mpi_kg_world, ierr )
!!$       deallocate(wk_zaj)
!!$
!!$       do nb = 0, nrank_e-1
!!$          do kb = 1, kimg
!!$             do jb = nis_e(nb), nie_e(nb)
!!$                do ib = 1, np_g1k(ik)
!!$                   zaj_ball(ib,jb,ik,kb) = wk_mpi(ib,jb-nis_e(nb)+1,kb,nb)
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$       deallocate(wk_mpi)
    end if
! ==============================================================================
#else
! ==============================================================================
    real(kind=DP) :: wk1(maxval(np_g1k),neg,kimg)
    integer :: ierr
    if(nrank_e == 1) then
       zaj_ball(:,:,ik,:) = zaj_l(:,:,ik,:)
    else
       wk1 = 0.0d0
       do i = 1, np_e
          wk1(:,neg_g(i),:) = zaj_l(:,i,ik,:)
       end do
       call mpi_allreduce(MPI_IN_PLACE,wk1,maxval(np_g1k)*neg*kimg &
            &            ,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       do i = 1, neg
          zaj_ball(:,i,ik,:) = wk1(:,i,:)
       end do
    end if
! ==============================================================================
#endif

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = 2
    end if
    allocate(p1Sp2(np_e,kimg_t))
    allocate(bp_w(nlmta,kimg_t))

    do i = 1, neg
       if(mod_pot == VDB) call broadcast_fs(i)  ! (fs(ri)_l(i,*,ik)-> bp_w
       call Psi1SPhi2_t(ik,i)                      ! -(contained here) -> p1Sp2
       call modify_bsd_and_phi_t(i)                ! -(c.h.) p1Sp2,phi_t -> bsd(ri),phi_t
    end do

    deallocate(p1Sp2)
    deallocate(bp_w)
  contains
#if 0
    subroutine broadcast_fs(io)
      integer, intent(in) :: io
      integer  :: ia, i, kimg_t

      i = map_z(io) ! MPI
      if(k_symmetry(ik) == GAMMA) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

      if(map_e(io) == myrank_e) then ! MPI
         if(k_symmetry(ik) == GAMMA) then
            do ia = 1, nlmta
               bp_w(ia,1) = fsr_l_2D(i,ia)
            end do
         else
            do ia = 1, nlmta
               bp_w(ia,1) = fsr_l_2D(i,ia)
               bp_w(ia,2) = fsi_l_2D(i,ia)
            end do
         end if
      end if ! MPI
      call mpi_bcast(bp_w,nlmta*kimg_t,mpi_double_precision,map_e(io),mpi_kg_world,ierr) ! MPI
    end subroutine broadcast_fs
#else
    subroutine broadcast_fs(io)
      integer, intent(in) :: io
      integer  :: ia, i, kimg_t
      integer :: sendcount, recvcounts(0:nrank_g-1), displs(0:nrank_g-1)
      real(kind=DP),allocatable,dimension(:,:) :: sendbuf
      sendcount = np_fs
      do i = 0, nrank_g - 1
         recvcounts(i) = nel_fs(i)
         displs(i)     = nis_fs(i) - 1
      enddo

      i = map_z(io) ! MPI
      if(k_symmetry(ik) == GAMMA) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(sendbuf(np_fs,kimg_t))

      if(map_e(io) == myrank_e) then ! MPI
         if(k_symmetry(ik) == GAMMA) then
            do ia = 1, np_fs
               sendbuf(ia,1) = fsr_l(i,ia,ik)
            end do
            call mpi_allgatherv(sendbuf(1,1),sendcount,        MPI_DOUBLE_PRECISION, &
                                bp_w(1,1),   recvcounts,displs,MPI_DOUBLE_PRECISION,mpi_ke_world,ierr)
         else
            do ia = 1, np_fs
               sendbuf(ia,1) = fsr_l(i,ia,ik)
               sendbuf(ia,2) = fsi_l(i,ia,ik)
            end do
            call mpi_allgatherv(sendbuf(1,1),sendcount,        MPI_DOUBLE_PRECISION, &
                                bp_w(1,1),   recvcounts,displs,MPI_DOUBLE_PRECISION,mpi_ke_world,ierr)
            call mpi_allgatherv(sendbuf(1,2),sendcount,        MPI_DOUBLE_PRECISION, &
                                bp_w(1,2),   recvcounts,displs,MPI_DOUBLE_PRECISION,mpi_ke_world,ierr)
         end if
      end if ! MPI
      call mpi_bcast(bp_w,nlmta*kimg_t,mpi_double_precision,map_e(io),mpi_kg_world,ierr) ! MPI
      deallocate(sendbuf)
    end subroutine broadcast_fs
#endif
    subroutine Psi1SPhi2_t(ik,i)
      integer, intent(in) :: ik,i

      real(kind=DP),allocatable,dimension(:,:)  :: p1Sp2_w

      integer   :: j,ia,p,q, kimg_t, iadd
      real(DP)  :: ar,ai

!!    DEBUG asms
!!    if((k_symmetry(1) == GAMMA .and. kimg == 2).or.kimg==1) then
      if((k_symmetry(ik) == GAMMA .and. kimg == 2).or.kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

      if(mod_pot == VANDERBILT_TYPE) then
         allocate(p1Sp2_w(np_e,kimg_t))
      end if

      p1Sp2 = 0.d0
      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, np_e ! MPI
            if(to == OTHER_BANDS .and. neg_g(j) == i) cycle ! MPI
            if(to == SAME_BAND .and. neg_g(j) /= i) cycle
            ar = 0.d0; if(kimg == 2) ai = 0.d0
            if(kimg == 1) then
               do ia = 1, nac
                  p = nlmta1(ia);      q = nlmta2(ia)
                  ar = ar + fqwei(ia)*(bp_w(p,1)*bsdr_l(j,q,ik) + bp_w(p,2)*bsdi_l(j,q,ik))
               end do
               p1Sp2(j,1) = ar
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, nac
                     p = nlmta1(ia);      q = nlmta2(ia)
                     ar = ar + fqwei(ia)*(bp_w(p,1)*bsdr_l(j,q,ik))
                  end do
                  p1Sp2(j,1) = ar
               else
                  do ia = 1, nac
                     p = nlmta1(ia);      q = nlmta2(ia)
                     ar = ar + fqwei(ia)*(bp_w(p,1)*bsdr_l(j,q,ik) + bp_w(p,2)*bsdi_l(j,q,ik))
                     ai = ai + fqwei(ia)*(bp_w(p,1)*bsdi_l(j,q,ik) - bp_w(p,2)*bsdr_l(j,q,ik))
                  end do
                  p1Sp2(j,1) = ar; p1Sp2(j,2) = ai
               end if
            end if
         end do
         p1Sp2_w = p1Sp2
         p1Sp2 = 0.d0
      end if

      do j = 1, np_e
         if(to == OTHER_BANDS .and. neg_g(j) == i) cycle
         if(to == SAME_BAND .and. neg_g(j) /= i) cycle
         if(kimg == 1) then
            do ia = ista_g1k(ik), iend_g1k(ik)
               iadd = ia - ista_g1k(ik) + 1
               p1Sp2(j,1) = p1Sp2(j,1) + zaj_ball(iadd,i,ik,1)*phi_l(iadd,j,ik,1)
            end do
            call mpi_allreduce(MPI_IN_PLACE,p1Sp2(j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = max(2,ista_g1k(ik)), iend_g1k(ik)
                  iadd = ia - ista_g1k(ik) + 1
                  ar = zaj_ball(iadd,i,ik,1)
                  ai = zaj_ball(iadd,i,ik,2)
                  p1Sp2(j,1) = p1Sp2(j,1) + (ar*phi_l(iadd,j,ik,1) + ai*phi_l(iadd,j,ik,2))*2.d0
               end do
               if(ista_g1k(ik) == 1) then
! ==== DEBUG by tkato 2011/09/03 ===============================================
!                 p1Sp2(j,1) = p1Sp2(j,1) + zaj_ball(1,i,ik,1)*zaj_ball(1,neg_g(j),ik,1)
                  p1Sp2(j,1) = p1Sp2(j,1) + zaj_ball(1,i,ik,1)*phi_l(1,j,ik,1)
! ==============================================================================
               endif
               call mpi_allreduce(MPI_IN_PLACE,p1Sp2(j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
            else
               do ia = ista_g1k(ik), iend_g1k(ik)
                  iadd = ia - ista_g1k(ik) + 1
                  ar = zaj_ball(iadd,i,ik,1)
                  ai = zaj_ball(iadd,i,ik,2)
                  p1Sp2(j,1) = p1Sp2(j,1) + ar*phi_l(iadd,j,ik,1) + ai*phi_l(iadd,j,ik,2)
                  p1Sp2(j,2) = p1Sp2(j,2) + ar*phi_l(iadd,j,ik,2) - ai*phi_l(iadd,j,ik,1)
               end do
               call mpi_allreduce(MPI_IN_PLACE,p1Sp2(j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
               call mpi_allreduce(MPI_IN_PLACE,p1Sp2(j,2),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
            end if
         end if
      end do
      if(mod_pot == VANDERBILT_TYPE) then                                ! MPI
         p1Sp2 = p1Sp2_w + p1Sp2                                         ! MPI
      endif

      if(mod_pot == VANDERBILT_TYPE) then
         deallocate(p1Sp2_w)
      end if

    end subroutine Psi1SPhi2_t

    subroutine modify_bsd_and_phi_t(i)
      integer, intent(in) :: i

      integer             :: j, ia, iadd
      real(DP)            :: sr, si

      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, np_e ! MPI
            if(to == OTHER_BANDS .and. neg_g(j) == i) cycle
            if(to == SAME_BAND .and. neg_g(j) /= i) cycle
            if(kimg == 1) then
               do ia = 1, nlmta
                  bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j,1)*bp_w(ia,1)
                  bsdi_l(j,ia,ik) = bsdi_l(j,ia,ik)-p1Sp2(j,1)*bp_w(ia,2)
               end do
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, nlmta
                     sr = bp_w(ia,1)
                     bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j,1)*sr
                  end do
               else
                  do ia = 1, nlmta
                     sr = bp_w(ia,1);  si = bp_w(ia,2)
                     bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                     bsdi_l(j,ia,ik) = bsdi_l(j,ia,ik)-p1Sp2(j,1)*si-p1Sp2(j,2)*sr
                  end do
               end if
            end if
         end do
      end if

      do j = 1, np_e
         if(to == OTHER_BANDS .and. neg_g(j) == i) cycle
         if(to == SAME_BAND .and. neg_g(j) /= i) cycle
         if(kimg == 1) then
            do ia = ista_g1k(ik), iend_g1k(ik)
               iadd = ia - ista_g1k(ik) + 1
               phi_l(iadd,j,ik,1) = phi_l(iadd,j,ik,1) - p1Sp2(j,1)*zaj_ball(iadd,i,ik,1)
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = ista_g1k(ik), iend_g1k(ik)
                  iadd = ia - ista_g1k(ik) + 1
                  sr = zaj_ball(iadd,i,ik,1) ;  si = zaj_ball(iadd,i,ik,2)
                  phi_l(iadd,j,ik,1) = phi_l(iadd,j,ik,1) - p1Sp2(j,1)*sr
                  phi_l(iadd,j,ik,2) = phi_l(iadd,j,ik,2) - p1Sp2(j,1)*si
               end do
            else
               do ia = ista_g1k(ik), iend_g1k(ik)
                  iadd = ia - ista_g1k(ik) + 1
                  sr = zaj_ball(iadd,i,ik,1) ;  si = zaj_ball(iadd,i,ik,2)
                  phi_l(iadd,j,ik,1) = phi_l(iadd,j,ik,1) - p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                  phi_l(iadd,j,ik,2) = phi_l(iadd,j,ik,2) - p1Sp2(j,1)*si-p1Sp2(j,2)*sr
               end do
            end if
         end if
      end do
    end subroutine modify_bsd_and_phi_t

  end subroutine mgs_sd2wf_each_k_G_3D


  subroutine m_ES_MGS_4_each_k(nfout,ik,mode)
    integer, intent(in) :: nfout,mode,ik ! mode={ORTHONORMALIZATION | NORMALIZATION}
    integer :: id_sname = -1,i
                                                  __TIMER_SUB_START(501)
#ifdef __FAPP__
    call fapp_start('mgs_4_each_k',1,1)
#endif
    call wd_title_of_the_operation()                        !-(c.h.)

!!$    call m_ES_mgs_alloc(ik)
    if(mode == ORTHONORMALIZATION) then
       call tstatc0_begin('modified_gram_schmidt ',id_sname,1)
    else
       call tstatc0_begin('modified_gram_schmidt (normalization only)',id_sname,1)
    endif
    if(modnrm == EXECUT) then
       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
       call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode,fsr_l,fsi_l,mod_pot=VDB)!-(m_E.S.)
       if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
       call tstatc0_end(id_sname)
!!$       if(mode==NORMALIZATION) call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
    else
       call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode,mod_pot=NORMCONSERVATION)!-(m_E.S.)
       if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
       call tstatc0_end(id_sname)
       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
    end if
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(:,ik) = OLD
#endif
!!$    call m_ES_mgs_dealloc
                                                  __TIMER_SUB_STOP(501)
#ifdef __FAPP__
    call fapp_stop('mgs_4_each_k',1,1)
#endif
  contains
    subroutine wd_title_of_the_operation
      integer, save :: iflag = 0
       if(iflag == 0 .and. ipri >= 2 .and. ik == 1) then
          if(modnrm == EXECUT) then
             write(nfout,*) ' <<< modified_gram_schmidt_vanderbilt_type >>>'
          else
             write(nfout,*) ' <<< modified_gram_schmidt_norm_conserve >>>'
          end if
       end if
       if(iflag == 0) iflag = 1
     end subroutine wd_title_of_the_operation
  end subroutine m_ES_MGS_4_each_k


!!$!BRANCH_P 3D_Parallel
!!$  subroutine m_ES_MGS_4_each_k_3D(nfout,ik,mode)
!!$    integer, intent(in) :: nfout,mode,ik ! mode={ORTHONORMALIZATION | NORMALIZATION}
!!$    integer :: id_sname = -1,i
!!$#ifdef __TIMER_SUB__
!!$  call timer_sta(502)
!!$#endif
!!$
!!$    call wd_title_of_the_operation()                        !-(c.h.)
!!$
!!$       call tstatc0_begin('modified_gram_schmidt ',id_sname,1)
!!$    if(modnrm == EXECUT) then
!!$       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
!!$       call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode &
!!$      &                      ,fsr_l,fsi_l,mod_pot=VDB)!-(m_E.S.)
!!$       call tstatc0_end(id_sname)
!!$    else
!!$       call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode &
!!$      &                      ,mod_pot=NORMCONSERVATION)!-(m_E.S.)
!!$       call tstatc0_end(id_sname)
!!$       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
!!$    end if
!!$
!!$#ifdef __TIMER_SUB__
!!$  call timer_end(502)
!!$#endif
!!$  contains
!!$    subroutine wd_title_of_the_operation
!!$      integer, save :: iflag = 0
!!$       if(iflag == 0 .and. ipri >= 1 .and. ik == 1) then
!!$          if(modnrm == EXECUT) then
!!$             write(nfout,*) ' <<< modified_gram_schmidt_vanderbilt_type >>>'
!!$          else
!!$             write(nfout,*) ' <<< modified_gram_schmidt_norm_conserve >>>'
!!$          end if
!!$       end if
!!$       if(iflag == 0) iflag = 1
!!$     end subroutine wd_title_of_the_operation
!!$  end subroutine m_ES_MGS_4_each_k_3D
!!$!BRANCH_P_END 3D_Parallel


  subroutine mgs_4_each_k_G_3D(k1,k2,ik,psi_l,mode,bpr_l,bpi_l,mod_pot,dryrun)
! Revised by T. Yamasaki in April 2006

    integer, intent(in)                                :: k1,k2,ik
!-F    real(kind=DP),dimension(kg1,np_e,k1:k2,kimg)      :: psi_l
    real(kind=DP),dimension(maxval(np_g1k),np_e,k1:k2,kimg)      :: psi_l
    integer, intent(in)                                :: mode,mod_pot
!-F    real(kind=DP),optional,dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l
    real(kind=DP),optional,dimension(np_e,np_fs,k1:k2) :: bpr_l,bpi_l
!-F PARA3D
    logical, intent(in), optional                      :: dryrun
    integer :: nbs, local_block, L_NB_STA, L_NB_END, in, ri, iy, ix, dsize_psi, dsize_bpri
    integer :: icount, iblk, jblk, i_NB, j_NB ,iq
    real(kind=DP),allocatable,dimension(:,:,:) :: wk_psi
    real(kind=DP),allocatable,dimension(:,:,:) :: wk_bpri
    real(kind=DP),allocatable,dimension(:,:,:) :: psi_t_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpr_t_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpi_t_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpr_tw1_dia
!!$    real(kind=DP),allocatable,dimension(:,:) :: bpr_tw2_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpi_tw1_dia
!!$    real(kind=DP),allocatable,dimension(:,:) :: bpi_tw2_dia
    real(kind=DP), allocatable, dimension(:,:,:,:) :: p1Sp2_t2_NB, p1Sp2_t1_NB

    logical :: drun
    integer ::       i
    real(kind=DP) :: fr
    integer  :: kimg_t_wk
    integer :: id_sname = -1
#ifdef MGS_DGEMM
    integer ::       NB_END, NB_END2, i1, i2
!!$    real(kind=DP), allocatable, dimension(:,:) :: bpr_tw1_BLAS
    real(kind=DP), allocatable, dimension(:,:,:) :: p1Sp2_NB
    integer, save :: ibsize_print = OFF
#endif
#ifndef SX
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache
    ncache = (m_CtrlP_cachesize()*1024)*3/4
! NEC tune <-------------------------------
#endif
                                                  __TIMER_SUB_START(503)

    call tstatc0_begin('mgs_4_each_k_G_3D ',id_sname)
    drun = .false.
    if(present(dryrun)) drun = dryrun

#ifdef MGS_DGEMM
    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
    if(ipri >= 2) then
       if(ibsize_print == OFF) then
          if(nblocksize_mgs_is_given) then
             write(nfout,'(" ! nblocksize_mgs_is_given")')
          else
             write(nfout,'(" ! nblocksize_mgs_is_given is false")')
          end if
          write(nfout,'( "! NB(=nblocksize_mgs) (mgs_4_each_k_G) = ",i8)') NB
          ibsize_print = ON
       end if
    end if
#endif

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!-F    call m_ES_mgs_alloc(ik)
!xx   call m_ES_mgs_alloc_3D(ik)
!!$    call m_ESortho_mgs_alloc_3D(ik)
    call m_ESortho_mgs_alloc(ik)
!!$    kimg_t_wk = kimg
!!$    if(kimg==2 .and. k_symmetry(ik) == GAMMA) kimg_t_wk = 1
    kimg_t_wk = 2
    if(k_symmetry(ik) == GAMMA) kimg_t_wk = 1
#ifdef MGS_DGEMM
    allocate(p1Sp2_NB(NB,NB,kimg_t_wk))
    allocate( wk_psi(mp_g1k(ik),NB,kimg) ) ; wk_psi = 0.0d0
    allocate( psi_t_dia(np_g1k(ik),NB,kimg) )
#endif
!PARA3D
    dsize_psi = mp_g1k(ik)*NB*kimg
    if(mod_pot == VANDERBILT_TYPE) then
      if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
!!$        i = 3
        i = 2
      else
!!$        i = 6
        i = 4
      endif
      allocate( wk_bpri(i,max(mp_fs,nac_p),NB) )
      dsize_bpri = i*max(mp_fs,nac_p)*NB
    endif
    allocate( bpr_t_dia(np_fs,NB) )
    allocate( bpi_t_dia(np_fs,NB) )
    allocate( bpr_tw1_dia(nac_p,NB) )
    allocate( bpi_tw1_dia(nac_p,NB) )

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_r_3D(k1,k2,ik,bpr_l,bpr_t)              ! bpr_l -> bpr_t
       else
          call m_ES_F_transpose_r_3D(k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)  ! bp[ri]_l -> bp[ri]_t
       end if
    end if

    call m_ES_W_transpose_r_3D(k1,k2,ik,psi_l,psi_t)    !-(m_E.S.) psi_ l-> psi_t
    !    This is, actually, not transpose.

     icount = (neg/NB+1)/nrank_e+1
!!$     allocate( p1Sp2_t2_NB(NB,NB,kimg_t_wk,icount) )
     allocate( p1Sp2_t1_NB(NB,NB,kimg_t_wk,icount) )

#ifdef MGS_DGEMM
    do i = 1, neg,NB
       NB_END = i + NB -1
       if( NB_END > neg ) NB_END = neg

       nbs = (i-1)/NB+1
       if(ipri>=2) &
            & write(nfout,'(" i, nbs = ",2i8, " myrank_e, lrank(nbs) = ",2i8)') i, nbs, myrank_e, lrank(nbs)
       if( myrank_e == lrank(nbs)) then
          local_block = nbsn(nbs)
          L_NB_STA = nbsn_sta(local_block)
          L_NB_END = nbsn_end(local_block)
          if(ipri>=2) write(nfout,'(" L_NB_STA, L_NB_END = ",2i8)') L_NB_STA, L_NB_END
!diagonal
                                                  __TIMER_DO_START(540)
    do i1 = L_NB_STA, L_NB_END
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,i1,mod_pot,fr,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t) ! fr = 1/dsqrt(<Psi(i)|S|Psi(i)>)
          if(dabs(fr-1.d0) > DELTA) &
               & call normalize_bp_and_psi_t_g(ik,i1,fr,mod_pot &
               & ,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t)
          !   |Psi(i)> = |Psi(i)> * fr,  <beta|Psi(i)> = <beta|Psi(i)> * fr
          if(mod_pot == VANDERBILT_TYPE) call cp_bpr2bprtw(i1,L_NB_END)  ! bpr_t -> bpr_tw1, bpr_tw2
       end if
       if(mode /= NORMALIZATION) then
          if(i1 == neg) cycle
          call cp_psi2psii_g(ik,i1) ! psi_t(:,i1,:) -> psi_ir,psi_ii
          call W1SW2_t_r_g(ik,i1,L_NB_END,mod_pot,psi_t,np_g1k_x,np_e,kimg) ! -> p1Sp2
          if(mod_pot == VANDERBILT_TYPE) &
               & call cp_bpr2bpi_g(kimg_t_wk, i1,bpr_t,bpi_t)  ! -> bp_ir, bp_ii
          call modify_bp_and_psi_t_r_g(ik,i1,L_NB_END,mod_pot &
               & ,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t) ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do   ! i1-loop
                                                  __TIMER_DO_STOP(540)
! bcast
       call cp_psi_bpri2dias_g()

       endif ! myrank_e == lrank(nbs)

    if(nrank_e > 1 ) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,529)
       call mpi_bcast(wk_psi,dsize_psi,mpi_double_precision,lrank(nbs),mpi_kg_world,ierr)
                                                 __TIMER_COMM_STOP(529)
    endif
                                                 __TIMER_COMM_START(530)
       do ri = 1, kimg
       do iy = 1, NB
         do ix = 1, np_g1k(ik)
           psi_t_dia(ix,iy,ri) = wk_psi(ix,iy,ri)
       enddo; enddo;  enddo;
                                                 __TIMER_COMM_STOP(530)
                                                 __TIMER_DO_START(541)
       do ri = 1, kimg
!       do iy = i, NB_END
       do iy = nbs_sta(nbs), nbs_end(nbs)
         do ix = 1, np_g1k(ik)
           zaj_ball(ix,iy,ik,ri) = wk_psi(ix,iy-nbs_sta(nbs)+1,ri)
       enddo; enddo;  enddo;
                                                 __TIMER_DO_STOP(541)
       if(mod_pot == VANDERBILT_TYPE) then
          if(nrank_e > 1 ) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,531)
              call mpi_bcast(wk_bpri,dsize_bpri,mpi_double_precision,lrank(nbs),mpi_kg_world,ierr)
                                                 __TIMER_COMM_STOP(531)
          endif
          if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
                                                 __TIMER_COMM_START(532)
             do iy = 1, NB
                do ix = 1, np_fs
                   bpr_t_dia(ix,iy) = wk_bpri(1,ix,iy)
                enddo
                do ix = 1, nac_p
                   bpr_tw1_dia(ix,iy) = wk_bpri(2,ix,iy)
!!$                   bpr_tw2_dia(ix,iy) = wk_bpri(3,ix,iy)
                enddo
             enddo
                                                 __TIMER_COMM_STOP(532)
                                                 __TIMER_DO_START(542)
              do iy = nbs_sta(nbs), nbs_end(nbs)     !       do iy = i,NB_END
                do ix = 1,np_fs
                   fsr_ball(iy,ix,ik) = wk_bpri(1,ix,iy-nbs_sta(nbs)+1)
                enddo
             enddo
                                                 __TIMER_DO_STOP(542)
          else
                                                 __TIMER_COMM_START(533)
             do iy = 1, NB
                do ix = 1, np_fs
                   bpr_t_dia(ix,iy) = wk_bpri(1,ix,iy)
                   bpi_t_dia(ix,iy) = wk_bpri(2,ix,iy)
                enddo
                do ix = 1, nac_p
                   bpr_tw1_dia(ix,iy) = wk_bpri(3,ix,iy)
!!$                   bpr_tw2_dia(ix,iy) = wk_bpri(4,ix,iy)
                   bpi_tw1_dia(ix,iy) = wk_bpri(4,ix,iy)
!!$                   bpi_tw2_dia(ix,iy) = wk_bpri(6,ix,iy)
                enddo
             enddo
                                                 __TIMER_COMM_STOP(533)
                                                 __TIMER_DO_START(543)
             do iy = nbs_sta(nbs), nbs_end(nbs)      !      do iy = i,NB_END
                do ix = 1,np_fs
                   fsr_ball(iy,ix,ik) = wk_bpri(1,ix,iy-nbs_sta(nbs)+1)
                   fsi_ball(iy,ix,ik) = wk_bpri(2,ix,iy-nbs_sta(nbs)+1)
                enddo
             enddo
                                                 __TIMER_DO_STOP(543)
          endif
       endif

!lower1
    if(mode /= NORMALIZATION) then
    icount = 0
    p1Sp2_t1_NB = 0.0d0
                                                 __TIMER_DO_START(544)
    do i2 = i+NB, neg,NB
       if(NB_END == neg) cycle
       NB_END2 = i2 + NB -1
       if( NB_END2 > neg ) NB_END2 = neg

       nbs = (i2-1)/NB+1
       if( myrank_e /= lrank(nbs)) cycle
       local_block = nbsn(nbs)
       L_NB_STA = nbsn_sta(local_block)
       L_NB_END = nbsn_end(local_block)
       if(mod_pot == VANDERBILT_TYPE) call cp_bpr2bprtw(L_NB_STA,L_NB_END)  ! bpr_t -> bpr_tw1, bpr_tw2

       call W1SW2_t_r_block_g(ik,i,L_NB_STA,p1Sp2_NB,NB_END,L_NB_END,kimg_t_wk &
            & , mod_pot,psi_t,psi_t_dia,np_g1k_x,np_e,kimg,bpr_tw1_dia,bpi_tw1_dia) ! -> p1Sp2
       icount = icount + 1

!-$ tune FUJITSU for block ----->>
!modify2010
                                                 __TIMER_COMM_START(534)
        do i_NB = 1, NB
           do iq = 1, kimg_t_wk
             do j_NB = 1,NB
               p1Sp2_t1_NB(j_NB,i_NB,iq,icount) = p1Sp2_NB(j_NB,i_NB,iq)
         end do ; enddo ; enddo
                                                 __TIMER_COMM_STOP(534)
!      if(nrank_g > 1 ) then
!       do i_NB = 1, NB
!          do iq = 1, kimg_t_wk
!            do j_NB = 1,NB
!              p1Sp2_t1_NB(j_NB,i_NB,iq,icount) = p1Sp2_NB(j_NB,i_NB,iq)
!        end do ; enddo ; enddo
!      endif
!-$ tune FUJITSU for block <<-----

    end do   ! i2-loop1
                                                 __TIMER_DO_STOP(544)

    if (NB_END /= neg ) then

! allreduce
    if(nrank_g > 1 ) then
      ix =  NB*NB*kimg_t_wk*icount
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,535)
      call mpi_allreduce(MPI_IN_PLACE, p1Sp2_t1_NB,ix,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(535)
    endif

    endif   !(NB_END /= neg )
!lower2
    icount = 0
                                                 __TIMER_DO_START(545)
    do i2 = i+NB, neg,NB
       if(NB_END == neg) cycle
       NB_END2 = i2 + NB -1
       if( NB_END2 > neg ) NB_END2 = neg
       nbs = (i2-1)/NB+1
       if( myrank_e /= lrank(nbs)) cycle
       local_block = nbsn(nbs)
       L_NB_STA = nbsn_sta(local_block)
       L_NB_END = nbsn_end(local_block)

!modify2010
       if(nrank_g > 1 ) then
         icount =  icount+1
                                                 __TIMER_DO_START(546)
         do i_NB = 1, NB
           do iq = 1, kimg_t_wk
             do j_NB = 1, NB
               p1Sp2_NB(j_NB,i_NB,iq) = p1Sp2_t1_NB(j_NB,i_NB,iq,icount)
         end do ; enddo ; enddo
                                                 __TIMER_DO_STOP(546)
       else
         icount =  icount+1
                                                 __TIMER_DO_START(547)
         do i_NB = 1, NB
           do iq = 1, kimg_t_wk
             do j_NB = 1, NB
               p1Sp2_NB(j_NB,i_NB,iq) = p1Sp2_t1_NB(j_NB,i_NB,iq,icount)
         end do ; enddo ; enddo
                                                 __TIMER_DO_STOP(547)
       endif
!      if(nrank_g > 1 ) then
!        icount =  icount+1
!        do i_NB = 1, NB
!          do iq = 1, kimg_t_wk
!            do j_NB = 1, NB
!              p1Sp2_NB(j_NB,i_NB,iq) = p1Sp2_t2_NB(j_NB,i_NB,iq,icount)
!        end do ; enddo ; enddo
!      endif

       call modify_bp_and_psi_t_r_blk_g(ik,i,L_NB_STA,p1Sp2_NB,NB_END,L_NB_END,kimg_t_wk &
            & , mod_pot,psi_t,psi_t_dia,np_g1k_x,np_e,kimg,bpr_t,bpi_t,bpr_t_dia,bpi_t_dia)
       !                   psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
    end do   ! i2-loop2
                                                 __TIMER_DO_STOP(545)
    end if   ! lower1 (mode /= NORMALIZATION)
    end do   ! i-loop
#else
    do i = 1, neg
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,i,mod_pot,fr,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t) ! fr = 1/dsqrt(<Psi(i)|S|Psi(i)>)
          if(dabs(fr-1.d0) > DELTA) &
               & call normalize_bp_and_psi_t_g(ik,i,fr,mod_pot,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t)
          !   |Psi(i)> = |Psi(i)> * fr,  <beta|Psi(i)> = <beta|Psi(i)> * fr
       end if
       if(mode /= NORMALIZATION) then
          if(i == neg) cycle
          call cp_psi2psii_g(ik,i) ! psi_t(:,i,:) -> psi_ir,psi_ii
          if(mod_pot == VANDERBILT_TYPE) &
               & call cp_bpr2bpi_g(kimg_t_wk,i,bpr_t,bpi_t) ! -> bp_ir, bp_ii
          call W1SW2_t_r_g(ik,i,neg,mod_pot,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t) ! -> p1Sp2
          call modify_bp_and_psi_t_r_g(ik,i,neg,mod_pot &
               & ,psi_t,np_g1k_x,np_e,kimg,bpr_t,bpi_t) ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do
#endif

    call m_ES_W_transpose_back_r_3D(k1,k2,ik,psi_l,psi_t)
    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_back_r_3D(k1,k2,ik,bpr_l,bpr_t)
       else
          call m_ES_F_transpose_back_r_3D(k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
       end if
    end if
!!$    deallocate(bpi_tw2_dia,bpi_tw1_dia,bpr_tw2_dia,bpr_tw1_dia)
    deallocate(p1Sp2_t1_NB)
    deallocate(bpi_tw1_dia,bpr_tw1_dia)
    deallocate(bpi_t_dia,bpr_t_dia)
    if(mod_pot == VANDERBILT_TYPE) deallocate(wk_bpri)
#ifdef MGS_DGEMM
    deallocate(psi_t_dia)
    deallocate(wk_psi)
    deallocate(p1Sp2_NB)
#endif
#else
    integer ::       nmax, ito
    integer, allocatable, dimension(:)         :: ib2to_a, ib2back_a
    allocate(ib2to_a(neg)); allocate(ib2back_a(neg))

!xx call m_ES_mgs_alloc(ik)
!!$    call m_ESortho_mgs_alloc_3D(ik)
    call m_ESortho_mgs_alloc(ik)

!!$    if(mod_pot == VANDERBILT_TYPE) call m_ES_F_transpose(k1,k2,ik,bpr_l,bpi_l,bpr_t,bpi_t)
    if(mod_pot == VANDERBILT_TYPE) &
         & call m_ES_F_transpose_r(.false.,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
    call m_ES_W_transpose_r(.false.,k1,k2,ik,psi_l,psi_t)    !-(m_E.S.) psi_ l-> psi_t
    do i = 1, neg
       ito = neordr(i,ik)
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,ito,mod_pot,fr,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t) ! fr = 1/dsqrt(<Psi(ito)|S|Psi(ito)>)
          if(dabs(fr-1.d0) > DELTA) &
               & call normalize_bp_and_psi_t_g(ik,ito,fr,mod_pot &
               & ,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t) !
          !   |Psi(ito)> = |Psi(ito)> * fr,  <beta|Psi(ito)> = <beta|Psi(ito)> * fr
       end if
       if(mode /= NORMALIZATION) then
          call substitute_jto_ib2back(i,nmax) !-(c.h.) ->ib2to_a,ib2back_a,
          if(nmax == 0) cycle
          call W1SW2_t(i,ito) ! -> p1Sp2
          call modify_bp_and_psi_t(i,ito) ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do

    call m_ES_W_transpose_back_r(.false.,k1,k2,ik,psi_l,psi_t)
!!$    if(mod_pot == VANDERBILT_TYPE) call m_ES_F_transpose_back(k1,k2,ik,bpr_l,bpi_l,bpr_t,bpi_t)
    if(mod_pot == VANDERBILT_TYPE) &
         & call m_ES_F_transpose_back_r(.false.,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
    deallocate(ib2back_a,ib2to_a)
#endif
!xx call m_ES_mgs_dealloc_3D()
    call m_ESortho_mgs_dealloc()
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(503)

!!$#else
!!$    integer, parameter :: npfsr_o = 0
!!$    integer            :: npfsi_o, npzr_o, npzi_o
!!$    integer, dimension(2) :: npzri_o
!!$    real(kind=DP)      :: fr,fi
!!$    integer            :: j, ibaik
!!$    integer                                         :: nsize_of_brd_w
!!$    real(kind=DP),allocatable, dimension(:)         :: brd_w
!!$
!!$    call alloc_and_brd_w()  ! nsize_of_brd_w, brd_w
!!$
!!$    ibaik = iba(ik)
!!$
!!$    npfsi_o = nlmta
!!$    call set_npzri(mod_pot,ik,npzr_o,npzi_o,npzri_o)  ! ->(npzr_o,npzi_o,npzri_o)
!!$
!!$    if(ipri >= 2) print '(" <<< modified_gram_schmidt_vanderbilt_type(cyclic)>>>")'
!!$    do i = 1, neg
!!$       ito = neordr(i,ik)
!!$       if(map_e(ito) == myrank_e) then               ! MPI
!!$          if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
!!$             call WSW(ito,fr) !   fr = 1/dsqrt(<Psi(ito)|S|Psi(ito)>)
!!$             if(dabs(fr-1.d0) > DELTA) call normalize_bp_and_psi(ito,fr)
!!$             !   |Psi(ito)> = |Psi(ito)> * fr,  <beta|Psi(ito)> = <beta|Psi(ito)> * fr
!!$          end if
!!$          if(mode /= NORMALIZATION) &
!!$               & call cp_bp_and_psi_2_brd_w2(k1,k2,ik,ito,mod_pot &
!!$               & , bpr_l,bpi_l,psi_l,npfsr_o,npfsi_o,npzri_o,nsize_of_brd,brd_w)
!!$          !                         -(m_E.S.) (bpr_l,bpi_l),psi_l -> brd_w   MPI
!!$       end if
!!$       if(mode /= NORMALIZATION) then
!!$          call mpi_bcast(brd_w,nsize_of_brd_w,mpi_double_precision &
!!$               & ,map_e(ito),mpi_k_world(myrank_k),ierr) !MPI
!!$
!!$          do j = ista_e, iend_e, istep_e                ! MPI
!!$             if(nrvf_ordr(j,ik) <= i) cycle
!!$             call W1SW2(j,fr,fi)
!!$             call modify_bp_and_psi(j,fr,fi)
!!$          end do
!!$       end if
!!$    end do
!!$    call dealloc_and_brd_w
!!$#endif
  contains

#ifdef MGS_DGEMM
    subroutine cp_psi_bpri2dias_g()
      integer :: ri,i1,ix,iy
                                                 __TIMER_COMM_START(527)
      do ri = 1, kimg
         do i1 = L_NB_STA, L_NB_END
            iy = i1-L_NB_STA+1
            do ix = 1, np_g1k(ik)
               wk_psi(ix,iy,ri) = psi_t(ix,i1,ri)
            enddo
         enddo
      enddo
                                                 __TIMER_COMM_STOP(527)
                                                 __TIMER_COMM_START(528)
      if(mod_pot == VANDERBILT_TYPE) then
         if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
            do i1 = L_NB_STA, L_NB_END
               iy = i1-L_NB_STA+1
               do ix = 1, np_fs
                  wk_bpri(1,ix,iy) = bpr_t(ix,i1)
               enddo
               do ix = 1, nac_p
                  wk_bpri(2,ix,iy) = bpr_tw1(ix,i1)
!!$                     wk_bpri(3,ix,iy) = bpr_tw2(ix,i1)
               enddo
            enddo
         else
            do i1 = L_NB_STA, L_NB_END
               iy = i1-L_NB_STA+1
               do ix = 1, np_fs
                  wk_bpri(1,ix,iy) = bpr_t(ix,i1)
                  wk_bpri(2,ix,iy) = bpi_t(ix,i1)
               enddo
               do ix = 1, nac_p
                  wk_bpri(3,ix,iy) = bpr_tw1(ix,i1)
!!$                     wk_bpri(4,ix,iy) = bpr_tw2(ix,i1)
!!$                     wk_bpri(5,ix,iy) = bpi_tw1(ix,i1)
                  wk_bpri(4,ix,iy) = bpi_tw1(ix,i1)
!!$                     wk_bpri(6,ix,iy) = bpi_tw2(ix,i1)
               enddo
            enddo
         endif
      endif
                                                 __TIMER_COMM_STOP(528)
  end subroutine cp_psi_bpri2dias_g
#endif

#ifdef MGS_DGEMM
    subroutine cp_bpr2bprtw(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: i, ia, p, q, i_last
                                                  __TIMER_SUB_START(505)
!!$      if(i1==1) then
!!$         i_last = i2
!!$      else
!!$         i_last = i1
!!$      end if
#ifdef MGS_DGEMM_DEBUG
      i_last = i1
#else
      i_last = i2
#endif
      if(k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DO_START(552)
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = bpr_t(q,i)
            end do
         end do
                                                  __TIMER_DO_STOP(552)
      else
                                                  __TIMER_DO_START(553)
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = bpr_t(q,i)
               bpi_tw1(ia,i) = bpi_t(p,i)
               bpi_tw2(ia,i) = bpi_t(q,i)
            end do
         end do
                                                  __TIMER_DO_STOP(553)
      end if
                                                  __TIMER_SUB_STOP(505)
    end subroutine cp_bpr2bprtw
#endif

!!$    subroutine alloc_p1Sp2(n)
!!$      integer, intent(in) :: n
!!$      if(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2) then
!!$         kimg_t = 1
!!$      else
!!$         kimg_t = 2
!!$      end if
!!$      allocate(p1Sp2(n,kimg_t))
!!$    end subroutine alloc_p1Sp2
!!$
!!$    subroutine dealloc_p1Sp2
!!$      deallocate(p1Sp2)
!!$    end subroutine dealloc_p1Sp2

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT

#else

    subroutine substitute_jto_ib2back(i,nmax)
      integer, intent(in)  :: i
      integer, intent(out) :: nmax

      integer              ::jto, j
      jto = 0
      do j = 1, neg
         if(nrvf_ordr(j,ik) <= i) cycle
         jto = jto + 1
         ib2to_a(jto)  = j
         ib2back_a(j)  = jto
      end do
      nmax = jto
      if(ipri >= 2 .and. ik == 1) then
         write(nfout,'(" <<substitute_jto_ib2back>>")')
         do j = 1, neg
            write(nfout,'(" ib2back_a(",i3,") = ",i3)') j, ib2back_a(j)
         end do
         do j = 1, nmax
            write(nfout,'(" ib2to_a(",i3,") = ",i3)') j, ib2to_a(j)
         end do
      end if
    end subroutine substitute_jto_ib2back

    subroutine W1SW2_t(i,ito)
      integer, intent(in)        :: i,ito

      integer       :: j, ia, jto, p, q, kimg_t
      real(kind=DP) :: ar, ai
      real(kind=DP), allocatable, dimension(:,:)  :: p1Sp2_t1, p1Sp2_t2             ! MPI

!!$      call alloc_p1Sp2(nmax)
      if(npes > 1) then
         if((k_symmetry(ik) == GAMMA .and. kimg==2) .or. kimg == 1) then
            kimg_t = 1
         else
            kimg_t = 2
         end if
         allocate(p1Sp2_t2(nmax,kimg_t)); p1Sp2_t2 = 0.d0                  ! MPI
         allocate(p1Sp2_t1(nmax,kimg_t)); p1Sp2_t1 = 0.d0                  ! MPI
      end if

      p1Sp2 = 0.d0

      if(mod_pot == VANDERBILT_TYPE) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(jto,ia,ar,ai,p,q)
#endif
         do jto = 1, nmax                                                 ! MPI
            j = ib2to_a(jto)                                              ! MPI
            if(kimg == 1) then
               ar = 0.d0
               do ia = 1, nac_p
                  p = nlmta1_p(ia);            q = nlmta2_p(ia)
                  ar = ar + fqwei_p(ia)*(bpr_t(p,ito)*bpr_t(q,j)+bpi_t(p,ito)*bpi_t(q,j))
               end do
               p1Sp2(jto,1) = ar
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  ar = 0.d0
                  do ia = 1, nac_p
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     ar = ar + fqwei_p(ia)*(bpr_t(p,ito)*bpr_t(q,j))
                  end do
                  p1Sp2(jto,1) = ar
               else
                  ar = 0.d0; ai = 0.d0
                  do ia = 1, nac_p
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     ar = ar + fqwei_p(ia)*(bpr_t(p,ito)*bpr_t(q,j)+bpi_t(p,ito)*bpi_t(q,j))
                     ai = ai + fqwei_p(ia)*(bpr_t(p,ito)*bpi_t(q,j)-bpi_t(p,ito)*bpr_t(q,j))
                  end do
                  p1Sp2(jto,1) = ar;  p1Sp2(jto,2) = ai
               end if
            end if
         end do
      end if

      if(ipri >= 2 .and. ik == 1) then
         if(nmax > 1) then
            write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3, " ito = ", i3)') ik,i,ito
            write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,1),jto=1, nmax)
            if(kimg == 2 .and. k_symmetry(ik) /= GAMMA) &
                 & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,kimg),jto=1, nmax)
         end if
      end if

      if(kimg == 1) then
         do jto = 1, nmax
            j = ib2to_a(jto)
            do ia = 1, np_g1k(ik)                                           ! MPI
               p1Sp2(jto,1) = p1Sp2(jto,1) + psi_t(ia,ito,1)*psi_t(ia,j,1 ) ! MPI
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
            do jto = 1, nmax
               j = ib2to_a(jto)
               do ia = 2, np_g1k(ik)
                  ar  = psi_t(ia,ito,1)
                  ai  = psi_t(ia,ito,2)
                  p1Sp2(jto,1) = p1Sp2(jto,1)+(ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2))*2.d0 ! MPI
               end do
!!$               if(mype /= 0) then
               if(myrank_e /= 0) then
                  ar = psi_t(1,ito,1)
                  ai = psi_t(1,ito,2)
                  p1Sp2(jto,1) = p1Sp2(jto,1) + (ar*psi_t(1,j,1) + ai*psi_t(1,j,2))*2.d0
               else
                  p1Sp2(jto,1) = p1Sp2(jto,1) + psi_t(1,ito,1)*psi_t(1,j,1)
               end if
            end do
         else
            do jto = 1, nmax
               j = ib2to_a(jto)
               do ia = 1, np_g1k(ik)                                          ! MPI
                  ar  = psi_t(ia,ito,1)
                  ai  = psi_t(ia,ito,2)
                  p1Sp2(jto,1) = p1Sp2(jto,1)+ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2) ! MPI
                  p1Sp2(jto,2) = p1Sp2(jto,2)+ar*psi_t(ia,j,2)-ai*psi_t(ia,j,1) ! MPI
               end do
            end do
         end if
      end if
      if(npes > 1 ) then
         do q = 1, kimg_t
            do jto = 1, nmax
               p1Sp2_t1(jto,q) = p1Sp2(jto,q)
            end do
         end do
         call mpi_allreduce(p1Sp2_t1, p1Sp2_t2,nmax*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
         do q = 1, kimg_t
            do jto = 1, nmax
               p1Sp2(jto,q) = p1Sp2_t2(jto,q)
            end do
         end do
      end if

      if(ipri >= 2 .and. ik == 1) then
         if(nmax > 1) then
            write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3, " ito = ", i3)') ik,i,ito
            write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,1),jto=1, nmax)
            if(kimg == 2 .and. k_symmetry(ik) /= GAMMA) &
                 & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,kimg),jto=1, nmax)
         end if
      end if

      if(nrank_e*nrank_g > 1) deallocate(p1Sp2_t1,p1Sp2_t2)
    end subroutine W1SW2_t

    subroutine modify_bp_and_psi_t(i,ito)
      integer, intent(in) :: i, ito
      integer             :: j,ia,jto
      real(kind=DP) :: sr, si

      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, neg
            if(nrvf_ordr(j,ik) <= i) cycle
            jto = ib2back_a(j)
            if(kimg == 1) then
               do ia = 1, np_fs
                  bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(jto,1)*bpr_t(ia,ito)
                  bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(jto,1)*bpi_t(ia,ito)
               end do
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, np_fs
                     bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(jto,1)*bpr_t(ia,ito)
                  end do
               else
                  do ia = 1, np_fs
                     sr  =  bpr_t(ia,ito);     si  =  bpi_t(ia,ito)
                     bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(jto,1)*sr+p1Sp2(jto,2)*si
                     bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(jto,1)*si-p1Sp2(jto,2)*sr
                  end do
               end if
            end if
         end do
      end if

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!CDIR NOSYNC
!CDIR CONCUR(BY=1)
#endif
      do j = 1,neg
         if(nrvf_ordr(j,ik) <= i) cycle
         jto = ib2back_a(j)         ! <-- 99Jan19 T. Y.
         if(kimg == 1) then
            do ia = 1, np_g1k(ik)                 ! MPI
               psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(jto,1)*psi_t(ia,ito,1 )
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = 1, np_g1k(ik)
                  sr  =  psi_t(ia,ito,1 );   si  =  psi_t(ia,ito,2 )
                  psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(jto,1)*sr
                  psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(jto,1)*si
               end do
            else
               do ia = 1, np_g1k(ik)                ! MPI
                  sr  =  psi_t(ia,ito,1 );   si  =  psi_t(ia,ito,2 )
                  psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(jto,1)*sr+p1Sp2(jto,2)*si
                  psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(jto,1)*si-p1Sp2(jto,2)*sr
               end do
            end if
         end if
      end do
!!$      call dealloc_p1Sp2()
    end subroutine modify_bp_and_psi_t
#endif

  end subroutine mgs_4_each_k_G_3D

!!$  subroutine set_npzri(mod_pot,ik,npzr_o,npzi_o,npzri_o)
!!$    integer, intent(in)  :: mod_pot,ik
!!$    integer, intent(out) :: npzr_o, npzi_o, npzri_o(2)
!!$
!!$    if(mod_pot == NORMCONSERVATION) then
!!$       npzr_o = 0
!!$    else
!!$       if(k_symmetry(ik) == GAMMA) then
!!$          npzr_o = nlmta
!!$       else
!!$          npzr_o = nlmta * 2
!!$       end if
!!$    end if
!!$    npzi_o = npzr_o + iba(ik)
!!$    npzri_o(1) = npzr_o; npzri_o(2) = npzi_o
!!$  end subroutine set_npzri

#ifdef TRANSPOSE
!!$#endif !TRANSPOSE_ORIGINAL


#else
!! #ifdef TRANSPOSE else

#endif
!! #ifdef TRANSPOSE end

  subroutine m_ES_F_transpose_r_3D(k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)

    integer, intent(in)                                                :: k1,k2,ik
!-F    real(kind=DP), intent(in), dimension(np_e,nlmta,k1:k2)             :: bpr_l
!-F    real(kind=DP), intent(out),dimension(np_fs_x,neg)                  :: bpr_t
!-F    real(kind=DP), intent(in), optional,dimension(np_e,nlmta,k1:k2)    :: bpi_l
!-F    real(kind=DP), intent(out),optional,dimension(np_fs_x,neg)         :: bpi_t
    real(kind=DP), intent(in), dimension(np_e,np_fs,k1:k2)             :: bpr_l
    real(kind=DP), intent(out),dimension(np_fs_x,np_e)                  :: bpr_t
    real(kind=DP), intent(in), optional,dimension(np_e,np_fs,k1:k2)    :: bpi_l
    real(kind=DP), intent(out),optional,dimension(np_fs_x,np_e)         :: bpi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r

    integer :: datasize, nb_proc, j, ix,iy, kimg_t, ifrom, ito
    integer :: pe_s,pe_r
    integer :: id_sname = -1
#ifdef MGS_DGEMM_DEBUG
    integer :: ia,p,q
#endif
                                                  __TIMER_SUB_START(511)

    if(k_symmetry(ik) == GAMMA) then ! assuming kimg == 2 when k_symmetry(ik) == GAMMA.
       kimg_t = 1
    else
       kimg_t = 2
    end if
                                                  __TIMER_DO_START(573)
    if(npes == 1) then
       if(kimg_t == 1) then
          if(np_fs_x > nlmta .or. neg > np_e) bpr_t = 0.d0
          do iy = 1, neg
!-F             ifrom = neordr(iy,ik)
             do ix = 1, nlmta
!-F                bpr_t(ix,iy) = bpr_l(ifrom,ix,ik)
                bpr_t(ix,iy) = bpr_l(iy,ix,ik)
             end do
          end do
       else if(kimg_t == 2) then
          if(np_fs_x > nlmta .or. neg > np_e)then
             bpr_t = 0.d0; bpi_t = 0.d0
          end if


          do iy = 1, neg
!-F             ifrom = neordr(iy,ik)
             do ix = 1, nlmta
!-F                bpr_t(ix,iy) = bpr_l(ifrom,ix,ik)
!-F                bpi_t(ix,iy) = bpi_l(ifrom,ix,ik)
                bpr_t(ix,iy) = bpr_l(iy,ix,ik)
                bpi_t(ix,iy) = bpi_l(iy,ix,ik)
             end do
          end do
       end if

#ifdef MGS_DGEMM_DEBUG
       goto 9999
#else
       return
#endif
    end if
                                                  __TIMER_DO_STOP(573)

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_F_transpose_r ',id_sname)

!-F add set work array for PARA3D
    bpr_t = 0.d0
    if(kimg_t == 2) bpi_t = 0.d0
                                                  __TIMER_DO_START(574)
    if(kimg_t == 1) then
       do iy = 1, np_e
          do ix = 1, np_fs
             bpr_t(ix,iy) = bpr_l(iy,ix,ik)
          end do
       end do
    else if(kimg_t == 2) then
       do iy = 1, np_e
          do ix = 1, np_fs
             bpr_t(ix,iy) = bpr_l(iy,ix,ik)
             bpi_t(ix,iy) = bpi_l(iy,ix,ik)
          end do
       end do
    end if
                                                  __TIMER_DO_STOP(574)
                                                  __TIMER_SUB_STOP(511)
  end subroutine m_ES_F_transpose_r_3D

  subroutine m_ES_F_transpose_back_r_3D(k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)

    integer, intent(in)                                                :: k1,k2,ik
!-F    real(kind=DP), intent(out),dimension(np_e,nlmta,k1:k2)             :: bpr_l
!-F    real(kind=DP), intent(in) ,dimension(np_fs_x,neg)                  :: bpr_t
!-F    real(kind=DP), intent(out),optional,dimension(np_e,nlmta,k1:k2)    :: bpi_l
!-F    real(kind=DP), intent(in) ,optional,dimension(np_fs_x,neg)         :: bpi_t
    real(kind=DP), intent(out), dimension(np_e,np_fs,k1:k2)             :: bpr_l
    real(kind=DP), intent(in),dimension(np_fs_x,np_e)                  :: bpr_t
    real(kind=DP), intent(out), optional,dimension(np_e,np_fs,k1:k2)    :: bpi_l
    real(kind=DP), intent(in),optional,dimension(np_fs_x,np_e)         :: bpi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r

    integer :: datasize, nb_proc, j, ix,iy, kimg_t,ito,ifrom
    integer :: pe_s,pe_r
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(512)

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = 2
    end if
                                                  __TIMER_DO_START(575)
    if(npes == 1) then
       if(kimg_t == 1) then
          do iy = 1, neg
!-F             ito = neordr(iy,ik)
             do ix = 1, nlmta
!-F                bpr_l(ito,ix,ik) = bpr_t(ix,iy)
                bpr_l(iy,ix,ik) = bpr_t(ix,iy)
             end do
          end do
       else
          do iy = 1, neg
!-F             ito = neordr(iy,ik)
             do ix = 1, nlmta
!-F                bpr_l(ito,ix,ik) = bpr_t(ix,iy)
!-F                bpi_l(ito,ix,ik) = bpi_t(ix,iy)
                bpr_l(iy,ix,ik) = bpr_t(ix,iy)
                bpi_l(iy,ix,ik) = bpi_t(ix,iy)
             end do
          end do
       end if
       return
    end if
                                                  __TIMER_DO_STOP(575)

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_F_transpose_back_r_3D ',id_sname)

!-F add return work array for PARA3D
                                                  __TIMER_DO_START(576)
    if(kimg_t == 1) then
       do iy = 1, np_e
          do ix = 1, np_fs
             bpr_l(iy,ix,ik) = bpr_t(ix,iy)
          end do
       end do
    else if(kimg_t == 2) then
       do iy = 1, np_e
          do ix = 1, np_fs
             bpr_l(iy,ix,ik) = bpr_t(ix,iy)
             bpi_l(iy,ix,ik) = bpi_t(ix,iy)
          end do
       end do
    end if
                                                  __TIMER_DO_STOP(576)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(512)
  end subroutine m_ES_F_transpose_back_r_3D

  subroutine m_ES_W_transpose_r_3D(k1,k2,ik,psi_l,psi_t)

    integer, intent(in)                                       :: k1,k2,ik
!-F    real(kind=DP), intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l  ! MPI
!-F    real(kind=DP), intent(out),dimension(np_g1k_x,neg,kimg)   :: psi_t  ! MPI
    real(kind=DP), intent(in), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l  ! MPI
    real(kind=DP), intent(out),dimension(np_g1k_x,np_e,kimg)   :: psi_t  ! MPI

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r

    integer :: ri, datasize, nb_proc, j, ix,iy, ifrom, ito
    integer :: pe_s,pe_r
    integer :: mp_g1k_x
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(513)
                                                  __TIMER_DO_START(577)
    if(npes == 1) then
       if(np_g1k(ik) > kg1 .or. neg > np_e) psi_t = 0.d0
       if(k_symmetry(ik) == GAMMA) then ! assuming kimg == 2 when k_symmetry(ik) == GAMMA.
          do iy = 1, neg
!-F             ifrom = neordr(iy,ik)
             do ix = 2, iba(ik)
!-F                psi_t(ix,iy,1) = psi_l(ix,ifrom,ik,1)
!-F                psi_t(ix,iy,2) = psi_l(ix,ifrom,ik,2)
                psi_t(ix,iy,1) = psi_l(ix,iy,ik,1)
                psi_t(ix,iy,2) = psi_l(ix,iy,ik,2)
             end do
             psi_t(1,iy,1) = psi_l(1,iy,ik,1)
             psi_t(1,iy,2) = 0.d0
          end do
       else
          do ri = 1, kimg
             do iy = 1, neg
!-F                ifrom = neordr(iy,ik)
                do ix = 1, iba(ik)
!-F                   psi_t(ix,iy,ri) = psi_l(ix,ifrom,ik,ri)
                   psi_t(ix,iy,ri) = psi_l(ix,iy,ik,ri)
                end do
             end do
          end do
       end if
       return
    end if
                                                  __TIMER_DO_STOP(577)

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_r_3D ',id_sname)

!-F add set work array for PARA3D
!    psi_t = 0.d0
                                                  __TIMER_DO_START(578)
    do ri = 1, kimg
       do iy = 1, np_e
          do ix = 1, np_g1k(ik)
             psi_t(ix,iy,ri) = psi_l(ix,iy,ik,ri)
          end do
       end do
    end do
                                                  __TIMER_DO_STOP(578)
                                                  __TIMER_DO_START(579)
    if(k_symmetry(ik) == GAMMA .and. myrank_g == 0 .and. kimg == 2) then
       do iy = 1, np_e
          psi_t(1,iy,2) = 0.d0
       end do
    end if
                                                  __TIMER_DO_STOP(579)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(513)
  end subroutine m_ES_W_transpose_r_3D

  subroutine m_ES_W_transpose_back_r_3D(k1,k2,ik,psi_l,psi_t)

    integer, intent(in)                                       :: k1,k2,ik
!-F    real(kind=DP), intent(out),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
!-F    real(kind=DP), intent(in) ,dimension(np_g1k_x,neg,kimg) :: psi_t
    real(kind=DP), intent(inout), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(in),dimension(np_g1k_x,np_e,kimg)   :: psi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r

    integer :: ri, datasize, nb_proc, j, ix,iy, ifrom,ito
    integer :: pe_s,pe_r
    integer :: mp_g1k_x
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(514)
                                                  __TIMER_DO_START(580)
    if(npes == 1) then
       do ri = 1, kimg
          do iy = 1, neg
!-F             ito = neordr(iy,ik)
             do ix = 1, iba(ik)
!-F                psi_l(ix,ito,ik,ri) = psi_t(ix,iy,ri)
                psi_l(ix,iy,ik,ri) = psi_t(ix,iy,ri)
             end do
          end do
       end do
       return
    end if
                                                  __TIMER_DO_STOP(580)

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_back_r ',id_sname)

!-F add return work array for PARA3D
                                                  __TIMER_DO_START(581)
    do ri = 1, kimg
       do iy = 1, np_e
          do ix = 1, np_g1k(ik)
             psi_l(ix,iy,ik,ri) = psi_t(ix,iy,ri)
          end do
       end do
    end do
                                                  __TIMER_DO_STOP(581)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(514)
  end subroutine m_ES_W_transpose_back_r_3D

! -- GLOBAL --
  subroutine WSW_t_g(ik,j,mod_pot,fr,psi_t,LA,LB,LC,bpr_t,bpi_t)
    integer, intent(in)        :: ik,j,mod_pot
    real(kind=DP), intent(out) :: fr
    integer, intent(in)        :: LA,LB,LC
    real(kind=DP),intent(in)   :: psi_t(LA,LB,LC)
    real(kind=DP),intent(in),optional :: bpr_t(:,:),bpi_t(:,:)

    integer              :: ia,p,q, i, ig1
    real(kind=DP)        :: fr1
    integer :: id_sname = -1
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('WSW_t_g ',id_sname)

                                                  __TIMER_SUB_START(504)
    fr = 0.d0
    if(mod_pot == VANDERBILT_TYPE) then
       if(k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DO_START(548)
          do ia = 1, nac_p
             p = nlmta1_p(ia);     q = nlmta2_p(ia)
             fr = fr+fqwei_p(ia)*(bpr_t(p,j)*bpr_t(q,j))
          end do
                                                  __TIMER_DO_STOP(548)
       else
                                                  __TIMER_DO_START(549)
          do ia = 1, nac_p
             p = nlmta1_p(ia);     q = nlmta2_p(ia)
             fr = fr+fqwei_p(ia)*(bpr_t(p,j)*bpr_t(q,j)+bpi_t(p,j)*bpi_t(q,j))
          end do
                                                  __TIMER_DO_STOP(549)
       end if
    end if

    fr1 = 0.d0
    if(kimg == 1) then
                                                  __TIMER_DO_START(550)
       do i = 1, np_g1k(ik)                       ! MPI
          fr1 = fr1 + psi_t(i,j,1)*psi_t(i,j,1)
       end do
                                                  __TIMER_DO_STOP(550)
    else if(kimg == 2) then
       ig1 = 1; if(k_symmetry(ik) == GAMMA .and. myrank_g == 0) ig1 = 2
                                                  __TIMER_DO_START(551)
       do i = ig1, np_g1k(ik)
          fr1 = fr1 + psi_t(i,j,1)*psi_t(i,j,1) + psi_t(i,j,2)*psi_t(i,j,2)
       end do
                                                  __TIMER_DO_STOP(551)
       if(k_symmetry(ik) == GAMMA) fr1 = fr1*2.d0
       if(ig1 == 2) fr1 = fr1 + (psi_t(1,j,1)*psi_t(1,j,1))
    end if
    fr = fr+fr1
    if(nrank_g > 1) then

                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,536)
       call mpi_allreduce(MPI_IN_PLACE,fr,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                  __TIMER_COMM_STOP(536)
    end if
    fr = 1.d0/dsqrt(fr)
    if(ipri >= 2 .and. (ik == 1 .or. ik == 2)) then
       write(nfout,'(" ((WSW_t_g)) ik = ",i8," fr = ",f21.15)') ik, fr
    end if
                                                  __TIMER_SUB_STOP(504)
  if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine WSW_t_g

! -- GLOBAL --
  subroutine normalize_bp_and_psi_t_g(ik,ibo,fr,mod_pot,psi_t,LA,LB,LC &
       & ,bpr_t,bpi_t)
    integer, intent(in)        :: ik,ibo
    real(kind=DP),intent(in)   :: fr
    integer, intent(in)        :: mod_pot,LA,LB,LC
    real(kind=DP),intent(inout),dimension(LA,LB,LC) :: psi_t
    real(kind=DP),intent(inout),optional :: bpr_t(:,:),bpi_t(:,:)

    integer                    :: ia,ri
    integer :: id_sname = -1

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('WSW_t_g ',id_sname)
                                                  __TIMER_SUB_START(506)
                                                  __TIMER_DO_START(554)
    if(mod_pot == VANDERBILT_TYPE) then
       if(k_symmetry(ik) == GAMMA) then
          do ia = 1, np_fs
             bpr_t(ia,ibo) = fr*bpr_t(ia,ibo)
          end do
       else
          do ia = 1, np_fs
             bpr_t(ia,ibo) = fr*bpr_t(ia,ibo)
             bpi_t(ia,ibo) = fr*bpi_t(ia,ibo)
          end do
       end if
    end if
                                                  __TIMER_DO_STOP(554)
                                                  __TIMER_DO_START(555)
    do ri = 1, kimg
!-F         psi_t(1:np_g1k(ik),ibo,ri) = fr * psi_t(1:np_g1k(ik),ibo,ri)
       psi_t(1:np_g1k(ik),ibo,ri) = fr * psi_t(1:np_g1k(ik),ibo,ri)
    end do
                                                  __TIMER_DO_STOP(555)
    if(ipri >= 2 .and. (ik == 1 .or. ik == 2)) then
       write(nfout,'(" ((normalize_bp_and_psi_t_g)) ik = ",i8," ibo = ",i8," fr = ",f21.15)') ik,ibo,fr
       write(nfout,'(6d14.6)') (psi_t(ia,ibo,1),ia=1, 6)
       if(kimg == 2) write(nfout,'(3d20.10)') (psi_t(ia,ibo,kimg),ia=1, 6)
    end if
                                                  __TIMER_SUB_STOP(506)
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine normalize_bp_and_psi_t_g

! -- GLOBAL --
  subroutine cp_bpr2bpi_g(kimg_t,i,bpr_t,bpi_t)
    integer, intent(in) :: kimg_t,i
    real(kind=DP),intent(in)           :: bpr_t(:,:)
    real(kind=DP),intent(in),optional  :: bpi_t(:,:)
    integer ::  nel
    integer :: id_sname = -1
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('cp_bpr2bpi_g ',id_sname)
    nel = min(np_fs,np_g1k_x)
    if(kimg_t == 1) then
       bp_ir(1:nel) = bpr_t(1:nel,i)
    else
       bp_ir(1:nel) = bpr_t(1:nel,i)
       bp_ii(1:nel) = bpi_t(1:nel,i)
    end if
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine cp_bpr2bpi_g

! -- GLOBAL --
  subroutine W1SW2_t_r_g(ik,i,NB_END,mod_pot,phi_t,LA,LB,LC &
!!$#ifndef MGS_DGEMM
       & ,phifr_t,phifi_t &
!!$#endif
       & )
    integer, intent(in)        :: ik,i, NB_END,mod_pot
    integer, intent(in)        :: LA,LB,LC
    real(kind=DP),intent(in)          :: phi_t(LA,LB,LC)
!!$#ifndef MGS_DGEMM
    real(kind=DP),intent(in),optional :: phifr_t(:,:),phifi_t(:,:)
!!$#endif

! Coded by T. Yamasaki in April 2006
! Revised according to the RIKEN phase tuning project 2009, 13 Sep 2009
! This subroutine is moved out from subroutine mgs_4_each_k_G() in a contained state.

    integer       :: j, ia, jto, p, q,  kimg_t
    real(kind=DP) :: ar, ai
    real(kind=DP), allocatable, dimension(:,:)  :: p1Sp2_t2, p1Sp2_t1  ! MPI
    character*4 F_RSVTASK
    integer       :: nt, mpant, mmdnt, ipar, ist, ied, mm

    integer :: id_sname = -1, id_sname2 = -1
    integer :: myrank_g_common
#ifndef SX
    integer :: n_unroll, jmax, ia_start
    integer :: ibsize,ibl1,ibl2
    integer :: ncache
                                                  __TIMER_SUB_START(507)
    ncache = (m_CtrlP_cachesize()*1024)*3/4
#endif

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_g ',id_sname)

    myrank_g_common = myrank_g


    if(nrank_e*nrank_g > 1) then
       if((k_symmetry(ik) == GAMMA .and. kimg==2) .or. kimg == 1) then
          kimg_t = 1
       else
          kimg_t = 2
       end if
       p = NB_END-i
       allocate(p1Sp2_t2(p,kimg_t)); p1Sp2_t2 = 0.d0                  ! MPI
       allocate(p1Sp2_t1(p,kimg_t)); p1Sp2_t1 = 0.d0                  ! MPI
    end if

    p1Sp2 = 0.d0

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg == 1) then
#ifndef SX
! NEC tune ------------------------------->
          if(ncache.eq.0) then
             ibsize=nac_p
!f
             if(ibsize == 0) ibsize = 1
          else
             ibsize=ncache/(8*(NB_END*2+3))
          endif
                                                  __TIMER_DO_START(556)
          do ibl1=1,nac_p,ibsize
             ibl2=ibl1+ibsize-1
             if(ibl2.gt.nac_p) ibl2=nac_p
! NEC tune <-------------------------------
#endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
!OCL NOFLTLD
             do j = i+1, NB_END
#ifdef SX
! NEC tune ------------------------------->
                ar = 0.d0
                do ia = 1, nac_p
#ifdef MGS_DGEMM
                   ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
#else
                   p = nlmta1_p(ia);            q = nlmta2_p(ia)
                   ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j)+bp_ii(p)*phifi_t(q,j))
#endif
                end do
                p1Sp2(j,1) = ar
#else
! NEC tune <-------------------------------
                ar = p1Sp2(j,1)
                do ia = ibl1, ibl2
#ifdef MGS_DGEMM
                   ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
#else
                   p = nlmta1_p(ia);            q = nlmta2_p(ia)
                   ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j)+bp_ii(p)*phifi_t(q,j))
#endif
                end do
                p1Sp2(j,1) = ar
#endif
            end do
#ifndef SX
! NEC tune
         end do
                                                  __TIMER_DO_STOP(556)
#endif
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
!!!$#ifdef NEC_TUNE_SMP
!!!$!CDIR PARALLEL DO private(ia,ar,p,q)
!!!$#endif
!!!$#ifdef VPP
!!!$*vocl loop, unroll(4)
!!!$#endif
!!!$!CDIR OUTERUNROLL=4
!!!$               do j = i+1, NB_END                                                ! MPI
!!!$                  ar = 0.d0
!!!$                  do ia = 1, nac_p
!!!$                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
!!!$                     ar = ar + fqwei_p(ia)*(bpr_t(p,i)*bpr_t(q,j))
!!!$                  end do
!!!$                  p1Sp2(j,1) = ar
!!!$               end do
#ifdef NEC_TUNE_SMP
             call getenv('F_RSVTASK',F_RSVTASK)
             read (F_RSVTASK,'(i4)') nt
#else
             nt = 1
#endif
             mpant = (NB_END-i)/nt
             mmdnt = mod(NB_END-i,nt)

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(mm,ist,ied,ia,j,p,q)
#endif
! NEC tune
!!#ifdef NEC_TUNE_SMP
#ifdef SX
             do ipar = 1, min(nt,NB_END-i)                                            ! SX
                IF (IPAR.LE.MMDNT) THEN                                               ! SX
                   MM = MPANT+1                                                       ! SX
                ELSE                                                                  ! SX
                   MM = MPANT                                                         ! SX
                ENDIF                                                                 ! SX
                IST = (IPAR-1)*MPANT + MIN(MMDNT+1,IPAR) + i                          ! SX
                IED = IST + mm - 1                                                    ! SX
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                do ia = 1, nac_p                                                      ! SX
                   do j = IST, IED
		     ! SX
#ifdef MGS_DGEMM
                      p1Sp2(j,1) = p1Sp2(j,1) + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))   ! SX
#else
                      p = nlmta1_p(ia);         q = nlmta2_p(ia)                      ! SX
                      p1Sp2(j,1) = p1Sp2(j,1) + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j))   ! SX
#endif
                   end do                                                             ! SX
                end do                                                                ! SX
             end do                                                                   ! SX
#else
! NEC tune ------------------------------->
             if(ncache.eq.0) then
                ibsize=nac_p
!f
                if(ibsize == 0) ibsize = 1
             else
                ibsize=ncache/(8*(NB_END*1+2))
             endif
                                                  __TIMER_DO_START(557)
!OCL NOFLTLD
             do ibl1=1,nac_p,ibsize
                ibl2=ibl1+ibsize-1
                if(ibl2.gt.nac_p) ibl2=nac_p
                do j = i+1, NB_END
                   if(ibl1.eq.1)then
                      ar=0.0d0
                   else
                      ar=p1Sp2(j,1)
                   endif
                   do ia = ibl1, ibl2
#ifdef MGS_DGEMM
                      ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))
#else
                      p = nlmta1_p(ia);         q = nlmta2_p(ia)
                      ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j))
#endif
                   end do
                   p1Sp2(j,1)=ar
                end do
             end do
                                                  __TIMER_DO_STOP(557)

! NEC tune <-------------------------------
#endif
          else      ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifndef SX
! NEC tune ------------------------------->
             if(ncache.eq.0) then
                ibsize=nac_p
!f
                if(ibsize == 0) ibsize = 1
             else
                ibsize=ncache/(8*(NB_END*4+3))
             endif
#ifdef DEBUG_MGS
             if(ipri >= 2) then
                write(nfout,'(" ibsize = ",i8," <<W1SW2_t_r>>")') ibsize
             end if
#endif
                                                  __TIMER_DO_START(558)
             do ibl1=1,nac_p,ibsize
                ibl2=ibl1+ibsize-1
                if(ibl2.gt.nac_p) ibl2=nac_p
! NEC tune <-------------------------------
#endif

#ifdef NEC_TUNE_SMP
#ifdef MGS_DGEMM
!CDIR PARALLEL DO private(ia,ar,ai)
#else
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
                do j = i+1, NB_END

#ifdef SX
                   ar = 0.d0; ai = 0.d0
                   do ia = 1, nac_p
#else
                   if(ibl1.eq.1) then
                      ar = 0.d0; ai = 0.d0
                   else
                      ar = p1Sp2(j,1)
                      ai = p1Sp2(j,2)

                   end if
                   do ia = ibl1, ibl2
#endif
#ifdef MGS_DGEMM
                      ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
                      ai = ai + fqwei_p(ia)*(bpr_tw1(ia,i)*bpi_tw2(ia,j)-bpi_tw1(ia,i)*bpr_tw2(ia,j))
#else
                      p = nlmta1_p(ia);         q = nlmta2_p(ia)
                      ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j)+bp_ii(p)*phifi_t(q,j))
                      ai = ai + fqwei_p(ia)*(bp_ir(p)*phifi_t(q,j)-bp_ii(p)*phifr_t(q,j))
#endif
                   end do
                   p1Sp2(j,1) = ar;  p1Sp2(j,2) = ai
                end do
#ifndef SX
! NEC tune
             end do
#endif
                                                  __TIMER_DO_STOP(558)
          end if
       end if
    end if
! --- <Psi_i|Psi_j> ---
    if(kimg == 1) then
#ifndef SX
! NEC tune ------------------------------------------------------------->
       if(ncache.eq.0) then
          ibsize=np_g1k(ik)
       else
          ibsize=ncache/(8*(NB_END*1+1))
       endif
                                                  __TIMER_DO_START(559)
       do ibl1=1,np_g1k(ik),ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.np_g1k(ik)) ibl2=np_g1k(ik)
! NEC tune <-------------------------------------------------------------
#endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
          do j = i+1, NB_END
#ifdef SX
! NEC tune ------------------------------------------------------------->
             do ia = 1, np_g1k(ik)                                           ! MPI
                p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(ia)*phi_t(ia,j,1 ) ! MPI
#else
             ai = p1Sp2(j,1)
             do ia = ibl1, ibl2
                ai = ai + psi_ir(ia)*phi_t(ia,j,1 ) ! MPI
! NEC tune <-------------------------------------------------------------
#endif
             end do

#ifndef SX
! NEC tune
             p1Sp2(j,1)=ai
#endif

          end do

#ifndef SX
! NEC tune
       end do
#endif
                                                  __TIMER_DO_STOP(559)

    else if(kimg == 2) then
       if(k_symmetry(ik) == GAMMA) then
!!!$            if(mype /= 0) then
!!!$#ifdef NEC_TUNE_SMP
!!!$!CDIR PARALLEL DO private(ia,ar,ai,p,q)
!!!$#endif
!!!$#ifdef VPP
!!!$*vocl loop, unroll(4)
!!!$#endif
!!!$!CDIR OUTERUNROLL=4
!!!$               do j = i+1, NB_END
!!!$                  do ia = 1, np_g1k(ik)
!!!$                     ar  = psi_t(ia,i,1)
!!!$                     ai  = psi_t(ia,i,2)
!!!$                     p1Sp2(j,1) = p1Sp2(j,1)+(ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2))*2.d0 ! MPI
!!!$                  end do
!!!$               end do
#ifdef NEC_TUNE_SMP
          call getenv('F_RSVTASK',F_RSVTASK)
          read (F_RSVTASK,'(i4)') nt
#else
          nt = 1
#endif

! NEC tune
!!#ifdef NEC_TUNE_SMP
#ifdef SX
          if( (NB_END-i)/nt .lt. 256 ) then

             if(myrank_g_common /= 0) then

                do j = i+1, NB_END
                   do ia = 1, np_g1k(ik)
                      ar  = psi_ir(ia)
                      ai  = psi_ii(ia)
                      p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                   end do
                end do
             else if(myrank_g_common == 0) then

                do j = i+1, NB_END
                   do ia = 2, np_g1k(ik)
                      ar  = psi_ir(ia)
                      ai  = psi_ii(ia)
                      p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                   end do
                   p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*phi_t(1,j,1)
                end do
             end if

          else

             mpant = (NB_END-i)/nt
             mmdnt = mod(NB_END-i,nt)
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(mm,ist,ied,ia,j,ar,ai)
#endif
                                                  __TIMER_DO_START(560)
             do ipar = 1, min(nt,NB_END-i)
                IF (IPAR.LE.MMDNT) THEN
                   MM = MPANT+1
                ELSE
                   MM = MPANT
                ENDIF
                IST = (IPAR-1)*MPANT + MIN(MMDNT+1,IPAR) + i
                IED = IST + mm - 1
                if(myrank_g_common /= 0) then
                     !     write(*,*)'ipar, IST, IED=', ipar, IST, IED
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                   do ia = 1, np_g1k(ik)
                      do j = IST, IED
                         ar  = psi_ir(ia)
                         ai  = psi_ii(ia)
                         p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                      end do
                   end do
                else if(myrank_g_common == 0) then
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                   do ia = 2, np_g1k(ik)
                      do j = IST, IED
                         ar  = psi_ir(ia)
                         ai  = psi_ii(ia)
                         p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                      end do
                   enddo
                   do j = IST, IED
                      p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*phi_t(1,j,1)
                   end do
                end if
             end do
                                                  __TIMER_DO_STOP(560)
          end if
! NEC tune ------------------------------------>
#else
!!$            if(ncache.eq.0) then
!!$               ibsize=np_g1k(ik)
!!$            else
!!$               ibsize=ncache/(8*(NB_END*2+2))
!!$            endif

!!$            allocate(psi_ir(np_g1k(ik)))
!!$            allocate(psi_ii(np_g1k(ik)))
!!$            do ia = 1, np_g1k(ik)
!!$               psi_ir(ia) = psi_t(ia,i,1)
!!$               psi_ii(ia) = psi_t(ia,i,2)
!!$            end do
!!$
!!$            n_unroll = 4
!!$            jto = (NB_END-i)/n_unroll
!!$            jmax = i + n_unroll*jto

          if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r(core) ',id_sname2)
          if(myrank_g_common /= 0) then
             ia_start = 1
          else if(myrank_g == 0) then
             ia_start = 2
          end if

                                                  __TIMER_DO_START(561)
          do j = i+1, NB_END
             ar=p1Sp2(j,1)
             do ia = ia_start, np_g1k(ik)
                ar = ar+(psi_ir(ia)*phi_t(ia,j,1)+psi_ii(ia)*phi_t(ia,j,2))*2.d0
             end do
             if(myrank_g_common == 0) then
                ar = ar + psi_ir(1)*phi_t(1,j,1)
             endif
             p1Sp2(j,1)=ar
          end do
                                                  __TIMER_DO_STOP(561)

!!$            do j = i+1, jmax, n_unroll
!!$               do ia = ia_start, np_g1k(ik)
!!$                  ar = psi_ir(ia)*2.d0; ai = psi_ii(ia)*2.d0
!!$                  p1Sp2(j  ,1) = p1Sp2(j  ,1) + (ar*psi_t(ia,j  ,1)+ai*psi_t(ia,j  ,2))
!!$                  p1Sp2(j+1,1) = p1Sp2(j+1,1) + (ar*psi_t(ia,j+1,1)+ai*psi_t(ia,j+1,2))
!!$                  p1Sp2(j+2,1) = p1Sp2(j+2,1) + (ar*psi_t(ia,j+2,1)+ai*psi_t(ia,j+2,2))
!!$                  p1Sp2(j+3,1) = p1Sp2(j+3,1) + (ar*psi_t(ia,j+3,1)+ai*psi_t(ia,j+3,2))
!!$               end do
!!$            end do
!!$            do j = jmax+1, neg
!!$               do ia = ia_start, np_g1k(ik)
!!$                  ar = psi_ir(ia); ai = psi_ii(ia)
!!$                  p1Sp2(j  ,1) = p1Sp2(j  ,1) + (ar*psi_t(ia,j,  1)+ai*psi_t(ia,j,  2))*2.d0
!!$               end do
!!$            end do
!!$            if(mype == 0) then
!!$               do j = i+1, neg
!!$                  p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*psi_t(1,j, 1)
!!$               end do
!!$            end if

          if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
#endif
! NEC tune <------------------------------------
       else    ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifdef SX
! NEC tune ------------------------------------>
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
                                                  __TIMER_DO_START(562)
          do j = i+1, NB_END
             do ia = 1, np_g1k(ik)                                          ! MPI
                ar  = psi_ir(ia)
                ai  = psi_ii(ia)
                p1Sp2(j,1) = p1Sp2(j,1)+ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2) ! MPI
                p1Sp2(j,2) = p1Sp2(j,2)+ar*phi_t(ia,j,2)-ai*phi_t(ia,j,1) ! MPI
             end do
          end do
                                                  __TIMER_DO_STOP(562)
#else
!!$            if(ncache.eq.0) then
          ibsize=np_g1k(ik)
!!$            else
!!$!!$               ibsize=ncache/(8*(neg*2+2))
!!$               ibsize=ncache/(8*(NB_END*2+2))
!!$            endif
!!$            allocate(psi_ir(np_g1k(ik)))
!!$            allocate(psi_ii(np_g1k(ik)))
!!$            do ia = 1, np_g1k(ik)
!!$               psi_ir(ia) = psi_iri(ia,1)
!!$               psi_ii(ia) = psi_iri(ia,2)
!!$            end do

                                                  __TIMER_DO_START(563)
          do ibl1=1,np_g1k(ik),ibsize
             ibl2=ibl1+ibsize-1
             if(ibl2.gt.np_g1k(ik)) ibl2=np_g1k(ik)
             if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r(core) ',id_sname2)
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#define _DIRECTIVE_UNROLLING_
#endif
#endif
#ifdef _DIRECTIVE_UNROLLING_
             do j = i+1, NB_END
                ar=p1Sp2(j,1)
                ai=p1Sp2(j,2)
                do ia = ibl1, ibl2
                   ar = ar+psi_ir(ia)*phi_t(ia,j,1)+psi_ii(ia)*phi_t(ia,j,2)
                   ai = ai+psi_ir(ia)*phi_t(ia,j,2)-psi_ii(ia)*phi_t(ia,j,1)
                end do
                p1Sp2(j,1)=ar
                p1Sp2(j,2)=ai
             end do
#else
             n_unroll = 4
             jto = (NB_END-i)/n_unroll
             jmax = i + n_unroll*jto
             do j = i+1, jmax, n_unroll
                do ia = ibl1, ibl2
                   ar = psi_ir(ia); ai = psi_ii(ia)
                   p1Sp2(j,1)   = p1Sp2(j,1)  +ar*phi_t(ia,j,1)  +ai*phi_t(ia,j,2)
                   p1Sp2(j,2)   = p1Sp2(j,2)  +ar*phi_t(ia,j,2)  -ai*phi_t(ia,j,1)
                   p1Sp2(j+1,1) = p1Sp2(j+1,1)+ar*phi_t(ia,j+1,1)+ai*phi_t(ia,j+1,2)
                   p1Sp2(j+1,2) = p1Sp2(j+1,2)+ar*phi_t(ia,j+1,2)-ai*phi_t(ia,j+1,1)
                   p1Sp2(j+2,1) = p1Sp2(j+2,1)+ar*phi_t(ia,j+2,1)+ai*phi_t(ia,j+2,2)
                   p1Sp2(j+2,2) = p1Sp2(j+2,2)+ar*phi_t(ia,j+2,2)-ai*phi_t(ia,j+2,1)
                   p1Sp2(j+3,1) = p1Sp2(j+3,1)+ar*phi_t(ia,j+3,1)+ai*phi_t(ia,j+3,2)
                   p1Sp2(j+3,2) = p1Sp2(j+3,2)+ar*phi_t(ia,j+3,2)-ai*phi_t(ia,j+3,1)
                end do
             end do
             do j = jmax+1, NB_END
                do ia = 1, ibl2-ibl1+1
                   ar = psi_ir(ia); ai=psi_ii(ia)
                   p1Sp2(j,1)   = p1Sp2(j,1)   + ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2)
                   p1Sp2(j,2)   = p1Sp2(j,2)   + ar*phi_t(ia,j,2)-ai*phi_t(ia,j,1)
                end do
             end do
#endif
             if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
          end do
                                                  __TIMER_DO_STOP(563)
!!$          deallocate(psi_ii, psi_ir)
! NEC tune <------------------------------------
#endif
       end if
    end if

    if(nrank_g > 1 .and. NB_END > i ) then
                                                  __TIMER_COMM_START(537)
       do q = 1, kimg_t
          do j = 1, NB_END-i
             p = j+i
             p1Sp2_t1(j,q) = p1Sp2(p,q)
          end do
       end do
                                                  __TIMER_COMM_STOP(537)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,538)
!!$         p = (neg-i)*kimg_t
!!$         call mpi_allreduce(p1Sp2_t1, p1Sp2_t2,(neg-i)*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
       call mpi_allreduce(p1Sp2_t1, p1Sp2_t2,(NB_END-i)*kimg_t,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                  __TIMER_COMM_STOP(538)
                                                  __TIMER_COMM_START(539)
       do q = 1, kimg_t
          do j = 1, NB_END-i
             p = j+i
             p1Sp2(p,q) = p1Sp2_t2(j,q)
          end do
       end do
!!$         p1Sp2(:,1:kimg_t) = p1Sp2_t2(:,1:kimg_t)
                                                  __TIMER_COMM_STOP(539)
    end if

    if(ipri >= 2 .and. ik == 1) then
       write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3)') ik,i
       write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(j,1),j=i+1,NB_END)
! === DEBUG by tkato 2013/11/05 ================================================
       if(nrank_e*nrank_g > 1) then
! ==============================================================================
       if(kimg_t == 2) &
            & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(j,kimg_t),j=i+1,NB_END)
! === DEBUG by tkato 2013/11/05 ================================================
       end if
! ==============================================================================
    end if
    if(nrank_e*nrank_g > 1) deallocate(p1Sp2_t1,p1Sp2_t2)
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(507)
  end subroutine W1SW2_t_r_g

! -- GLOBAL --
!-$ tune FUJITSU for block ----->>
!     13th Sep 2009
#ifdef MGS_DGEMM
  subroutine W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
       &  ,   mod_pot,psi_t,psi_t_dia,LA,LB,LC,bpr_tw1_dia,bpi_tw1_dia)
    integer, intent(in)        :: ik,i, i2, NB_END, NB_END2, kimg_t_wk
    integer, intent(in)        :: mod_pot,LA,LB,LC
    real(kind=DP), intent(out), dimension(NB,NB,kimg_t_wk) :: p1Sp2_NB
    real(kind=DP), intent(in), dimension(LA,LB,LC)      :: psi_t
    real(kind=DP), intent(in), dimension(LA,NB,LC)      :: psi_t_dia
    real(kind=DP), intent(in), dimension(nac_p,*),optional :: bpr_tw1_dia,bpi_tw1_dia

    real(kind=DP), allocatable, dimension(:,:) :: bpr_tw1_BLAS, bpi_tw1_BLAS ! work

    integer       :: iblk, jblk, i_NB, j_NB
    integer       :: ii, ia,  M, N, LDA, LDB, p

    integer :: id_sname = -1, id_sname2 = -1
                                                  __TIMER_SUB_START(509)
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_block_g ',id_sname)

    p1Sp2_NB = 0.0d0
    M = NB_END2-i2+1; N = NB_END-i+1
    if(mod_pot == VANDERBILT_TYPE .and. nac_p > 0) then
       allocate(bpr_tw1_BLAS(nac_p, NB))
       bpr_tw1_BLAS = 0.0
                                                  __TIMER_DO_START(570)
!OCL NOFLTLD
       do iblk = i, NB_END
          ii = iblk - i + 1
          do ia = 1, nac_p
             bpr_tw1_BLAS(ia,ii) = fqwei_p(ia)*bpr_tw1_dia(ia,ii)
          end do
       end do
                                                  __TIMER_DO_STOP(570)
       if(k_symmetry(ik) /= GAMMA) then
          allocate(bpi_tw1_BLAS(nac_p,NB))
          bpi_tw1_BLAS = 0.d0
                                                  __TIMER_DO_START(571)
          do iblk = i, NB_END
             ii = iblk - i + 1
             do ia = 1, nac_p
                bpi_tw1_BLAS(ia,ii) = fqwei_p(ia)*bpi_tw1_dia(ia,ii)
             end do
          end do
                                                  __TIMER_DO_STOP(571)
       end if

       if(kimg==1) then
                                                  __TIMER_DGEMM_START(515)
          call DGEMM__('T','N',M,N,nac_p,1.0d0,bpr_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
          call DGEMM__('T','N',M,N,nac_p,1.0d0,bpi_tw2(1,i2),nac_p,bpi_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
                                                  __TIMER_DGEMM_STOP(515)
!!$                     ar = ar + fqwei_p(ia)*(bpr_t(p,i)*bpr_t(q,j)+bpi_t(p,i)*bpi_t(q,j))
       else if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DGEMM_START(516)
          call DGEMM__('T','N',M,N,nac_p,1.0d0,bpr_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
!!$                     ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))
                                                  __TIMER_DGEMM_STOP(516)
       else if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
                                                  __TIMER_DGEMM_START(517)
          call DGEMM__('T','N',M,N,nac_p, 1.0d0,bpr_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
          call DGEMM__('T','N',M,N,nac_p, 1.0d0,bpi_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,2),NB)
          call DGEMM__('T','N',M,N,nac_p, 1.0d0,bpi_tw2(1,i2),nac_p,bpi_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
          call DGEMM__('T','N',M,N,nac_p,-1.0d0,bpr_tw2(1,i2),nac_p,bpi_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,2),NB)
!!$                     ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
!!$                     ai = ai + fqwei_p(ia)*(bpr_tw1(ia,i)*bpi_tw2(ia,j)-bpi_tw1(ia,i)*bpr_tw2(ia,j))
                                                  __TIMER_DGEMM_STOP(517)
       end if

       if(k_symmetry(ik) /= GAMMA) deallocate(bpi_tw1_BLAS)
       deallocate(bpr_tw1_BLAS)
    end if
! --------------------------------------
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_block_g(core) ',id_sname2)
    LDA = np_g1k_x; LDB = np_g1k_x
    if(kimg == 1) then
                                                  __TIMER_DGEMM_START(518)
       call DGEMM__('T','N',M,N,np_g1k(ik),1.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,1), LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
                                                  __TIMER_DGEMM_STOP(518)
    else if(kimg == 2 .and. k_symmetry(ik) == GAMMA ) then
                                                  __TIMER_DGEMM_START(519)
       call DGEMM__('T','N',M,N,np_g1k(ik),2.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,1),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik),2.0d0,psi_t(1,i2,2),LDA,psi_t_dia(1,1,2),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
                                                  __TIMER_DGEMM_STOP(519)
       if (myrank_g == 0 ) then
                                                  __TIMER_DO_START(572)
          do iblk = i, NB_END                !for BLAS3
             i_NB = iblk - i + 1
             do jblk = i2, NB_END2
                j_NB = jblk - i2 + 1
!f                  p1Sp2_NB(j_NB,i_NB,1) = p1Sp2_NB(j_NB,i_NB,1) &
!f                       & - (psi_t(1,iblk,1)*psi_t(1,jblk,1)+psi_t(1,iblk,2)*psi_t(1,jblk,2)*2.d0)
                p1Sp2_NB(j_NB,i_NB,1) = p1Sp2_NB(j_NB,i_NB,1) &
                     & - (psi_t_dia(1,i_NB,1)*psi_t(1,jblk,1)+psi_t_dia(1,i_NB,2)*psi_t(1,jblk,2)*2.d0)
             end do
          end do
                                                  __TIMER_DO_STOP(572)
       end if
    else  ! kimg==2 .and. k_symmetry(ik) /= GAMMA
                                                  __TIMER_DGEMM_START(520)
       call DGEMM__('T','N',M,N,np_g1k(ik), 1.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,1),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik), 1.0d0,psi_t(1,i2,2),LDA,psi_t_dia(1,1,2),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik), 1.0d0,psi_t(1,i2,2),LDA,psi_t_dia(1,1,1),LDB,1.0d0,p1Sp2_NB(1,1,2),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik),-1.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,2),LDB,1.0d0,p1Sp2_NB(1,1,2),NB)
                                                  __TIMER_DGEMM_STOP(570)
    end if
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)

!-F    if(nrank_e > 1 ) then
!-F       if(ipri>=2) then
!-F          write(nfout,'(" NB_END, NB_END2, NB_END-i+1, NB_END2-i2+1, NB, kimg_t_wk = ",6i8)') &
!-F               & NB_END, NB_END2,  NB_END-i+1, NB_END2-i2+1, NB, kimg_t_wk
!-F       end if
!-F       call mpi_allreduce(MPI_IN_PLACE,p1Sp2_NB,NB*NB*kimg_t_wk,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
!-F    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(509)
  end subroutine W1SW2_t_r_block_g
#endif
!-$ tune FUJITSU for block <<-----

! -- GLOBAL --
#ifdef MGS_DGEMM
    subroutine modify_bp_and_psi_t_r_blk_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
	&   , mod_pot,psi_t,psi_t_dia,LA,LB,LC,bpr_t,bpi_t,bpr_t_dia,bpi_t_dia)
      integer, intent(in) :: ik,i,i2,mod_pot,LA,LB,LC
      real(kind=DP), intent(in),  dimension(NB,NB,kimg_t_wk) :: p1Sp2_NB
      real(kind=DP), intent(inout), dimension(LA,LB,LC) :: psi_t
      real(kind=DP), intent(in), dimension(LA,NB,LC)    :: psi_t_dia
      real(kind=DP), intent(inout), optional :: bpr_t(:,:),bpi_t(:,:)
      real(kind=DP), intent(in), optional :: bpr_t_dia(:,:),bpi_t_dia(:,:)

      integer, intent(in) :: NB_END, NB_END2, kimg_t_wk
      integer       ::  N, K, LDA, LDB, LDC
      integer :: id_sname = -1
                                                  __TIMER_SUB_START(510)
      if(i == neg) return

      if(sw_timing_2ndlevel == ON) call tstatc0_begin('modify_bp_and_psi_t_r_blk_g ',id_sname)

      N = NB_END2-i2+1
      K = NB_END-i+1

      if(mod_pot == VANDERBILT_TYPE .and. np_fs.gt.0) then
         LDA = np_fs_x; LDB = NB; LDC = np_fs_x
                                                  __TIMER_DGEMM_START(521)
         call DGEMM__('N','T',np_fs,N,K,-1.0d0,bpr_t_dia(1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.0d0,bpr_t(1,i2),LDC)
                                                  __TIMER_DGEMM_STOP(521)

         if(kimg == 1 .or. (kimg==2 .and. k_symmetry(ik) /= GAMMA)) then
                                                  __TIMER_DGEMM_START(522)
            call DGEMM__('N','T',np_fs,N,K,-1.0d0,bpi_t_dia(1,1),LDA,p1Sp2_NB(1,1,1),LDB, 1.0d0,bpi_t(1,i2),LDC)
                                                  __TIMER_DGEMM_STOP(522)
         end if

         if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
                                                  __TIMER_DGEMM_START(523)
            call DGEMM__('N','T',np_fs,N,K, 1.0d0,bpi_t_dia(1,1),LDA,p1Sp2_NB(1,1,2),LDB,1.0d0,bpr_t(1,i2),LDC)
            call DGEMM__('N','T',np_fs,N,K,-1.0d0,bpr_t_dia(1,1),LDA,p1Sp2_NB(1,1,2),LDB,1.0d0,bpi_t(1,i2),LDC)
                                                  __TIMER_DGEMM_STOP(523)
         end if
      end if

      LDA = np_g1k_x; LDB = NB; LDC = np_g1k_x
      if(kimg == 1) then
                                                  __TIMER_DGEMM_START(524)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,1),LDC)
                                                  __TIMER_DGEMM_STOP(524)
      else if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DGEMM_START(525)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,1),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,2),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,2),LDC)
                                                  __TIMER_DGEMM_STOP(525)
      else  ! kimg == 2 .and. k_symmetry(ik) /= GAMMA
                                                  __TIMER_DGEMM_START(526)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,1),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,2),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,2),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K, 1.d0,psi_t_dia(1,1,2),LDA,p1Sp2_NB(1,1,2),LDB,1.d0,psi_t(1,i2,1),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,2),LDB,1.d0,psi_t(1,i2,2),LDC)
                                                  __TIMER_DGEMM_STOP(526)
      end if
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(510)
    end subroutine modify_bp_and_psi_t_r_blk_g
!-$ tune FUJITSU for block <<-----
#endif

! -- GLOBAL --
  subroutine modify_bp_and_psi_t_r_g(ik,i,NB_END,mod_pot &
       & ,psi_t,LA,LB,LC,bpr_t,bpi_t)
! Coded by T. Yamasaki in April 2006
! Revised according to the RIKEN phase tuning project 2009, 13 Sep 2009
    integer, intent(in) :: ik,i, NB_END,mod_pot
    integer, intent(in) :: LA,LB,LC
    real(kind=DP),intent(inout) :: psi_t(LA,LB,LC)
    real(kind=DP),intent(inout),optional :: bpr_t(:,:),bpi_t(:,:)

    integer             :: j,ia
    real(kind=DP) :: sr, si
#ifndef SX
    real(kind=DP) :: ar, ai
#endif
#ifdef MGS_DGEMM
    integer :: p, q
#endif
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(508)
    if(i == neg) return
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('modify_bp_and_psi_t_r_g ',id_sname)

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(564)
          do j = i+1, NB_END
#ifdef SX
             do ia = 1, np_fs
                bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(j,1)*bp_ir(ia)
                bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(j,1)*bp_ii(ia)
             end do
#else
             ar=p1Sp2(j,1)
             do ia = 1, np_fs
                bpr_t(ia,j) = bpr_t(ia,j) - ar*bp_ir(ia)
                bpi_t(ia,j) = bpi_t(ia,j) - ar*bp_ii(ia)
             end do
#endif

          end do  ! j-loop
                                                  __TIMER_DO_STOP(564)
       else if(kimg == 2) then
          if(k_symmetry(ik) == GAMMA) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(565)
             do j = i+1, NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
                do ia = 1, np_fs
                   bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(j,1)*bp_ir(ia)
                end do
#else
                ar=p1Sp2(j,1)
                do ia = 1, np_fs
                   bpr_t(ia,j) = bpr_t(ia,j) - ar*bp_ir(ia)
                end do
! NEC tune <----------------------------------------------------------
#endif
             end do
                                                  __TIMER_DO_STOP(565)
          else  ! kimg==2, k_symmetry(ik) /= GAMMA
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(566)
             do j = i+1, NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
                do ia = 1, np_fs
                   sr  =  bp_ir(ia);     si  =  bp_ii(ia)
                   bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                   bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(j,1)*si-p1Sp2(j,2)*sr
                end do
#else
                ar=p1Sp2(j,1)
                ai=p1Sp2(j,2)
                do ia = 1, np_fs
                   sr  =  bp_ir(ia);     si  =  bp_ii(ia)
                   bpr_t(ia,j) = bpr_t(ia,j) - ar*sr+ai*si
                   bpi_t(ia,j) = bpi_t(ia,j) - ar*si-ai*sr
                end do
#endif
             end do
                                                  __TIMER_DO_STOP(566)
! NEC tune <----------------------------------------------------------
          end if
       end if
    end if

    if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(567)
       do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
          do ia = 1, np_g1k(ik)                 ! MPI
!!$             psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_t(ia,i,1 )
             psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_ir(ia )
          end do
#else
          ar=p1Sp2(j,1)
          do ia = 1, np_g1k(ik)
!!$             psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_t(ia,i,1 )
             psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_ir(ia)
          end do
! NEC tune <----------------------------------------------------------
#endif
       end do
                                                  __TIMER_DO_STOP(567)
    else if(kimg == 2) then
       if(k_symmetry(ik) == GAMMA) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(568)
          do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
             do ia = 1, np_g1k(ik)
!!$                psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_t(ia,i,1 )
!!$                psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(j,1)*psi_t(ia,i,2 )
                psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_ir(ia)
                psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(j,1)*psi_ii(ia)
             end do
#else
             ar=p1Sp2(j,1)
             do ia = 1, np_g1k(ik)
!!$                psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_t(ia,i,1 )
!!$                psi_t(ia,j,2) = psi_t(ia,j,2) - ar*psi_t(ia,i,2 )
                psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_ir(ia)
                psi_t(ia,j,2) = psi_t(ia,j,2) - ar*psi_ii(ia)
             end do
! NEC tune <----------------------------------------------------------
#endif
          end do
                                                  __TIMER_DO_STOP(568)
       else   ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(569)
          do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
             do ia = 1, np_g1k(ik)                ! MPI
!!$                sr  =  psi_t(ia,i,1 );   si  =  psi_t(ia,i,2 )
                sr  =  psi_ir(ia);   si  =  psi_ii(ia)
                psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(j,1)*si-p1Sp2(j,2)*sr
             end do
#else
             ar=p1Sp2(j,1)
             ai=p1Sp2(j,2)
             do ia = 1, np_g1k(ik)
!!$                sr  =  psi_t(ia,i,1 );   si  =  psi_t(ia,i,2 )
                sr  =  psi_ir(ia);   si  =  psi_ii(ia)
                psi_t(ia,j,1) = psi_t(ia,j,1) - ar*sr+ai*si
                psi_t(ia,j,2) = psi_t(ia,j,2) - ar*si-ai*sr
             end do
#endif
! NEC tune <----------------------------------------------------------
          end do
                                                  __TIMER_DO_STOP(569)
       end if
    end if
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(508)
  end subroutine modify_bp_and_psi_t_r_g


! -- GLOBAL --
  subroutine cp_psi2psii_g(ik,i)
    integer, intent(in) :: ik,i
    psi_ir(1:np_g1k(ik)) = psi_t(1:np_g1k(ik),i,1)
    if(kimg==2) psi_ii(1:np_g1k(ik)) = psi_t(1:np_g1k(ik),i,kimg)
  end subroutine cp_psi2psii_g


  subroutine mgs_phi2wf_each_k_G(ik,phi_l,mode,bsdr_l,bsdi_l,mod_pot)
    integer, intent(in)             :: ik
!!$    real(kind=DP), intent(inout)    :: phi_l(kg1,np_e,ik:ik,kimg)
    real(kind=DP),intent(inout),dimension(maxval(np_g1k),np_e,ik:ik,kimg) :: phi_l
    integer, intent(in)             :: mode,mod_pot
    real(kind=DP), optional,intent(inout), dimension(np_e,np_fs,ik:ik) :: bsdr_l,bsdi_l
!!$    real(kind=DP), optional,intent(inout), dimension(np_e,nlmta,ik:ik) :: bsdr_l,bsdi_l

    integer :: nbs, local_block, L_NB_STA, L_NB_END, in, ri, iy, ix, dsize_psi, dsize_bpri
    integer :: icount, iblk, jblk, i_NB, j_NB ,iq
    real(kind=DP),allocatable,dimension(:,:,:) :: wk_psi
    real(kind=DP),allocatable,dimension(:,:,:) :: wk_bpri
    real(kind=DP),allocatable,dimension(:,:,:) :: psi_t_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpr_t_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpi_t_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpr_tw1_dia
    real(kind=DP),allocatable,dimension(:,:) :: bpi_tw1_dia
    real(kind=DP), allocatable, dimension(:,:,:,:) :: p1Sp2_t2_NB, p1Sp2_t1_NB
!!$    real(kind=DP), allocatable, dimension(:,:,:) :: bsdr_lc, bsdi_lc

    integer ::       i, j

    real(kind=DP) :: fr
    integer :: kimg_t_wk

    real(kind=DP),allocatable,dimension(:,:)    :: phifr_t, phifi_t
    real(kind=DP),allocatable, dimension(:,:,:) :: phi_t ! d(maxval(np_g1k),np_e,kimg)
#ifdef MGS_DGEMM
    integer ::       NB_END, NB_END2, i1, i2
    real(kind=DP), allocatable, dimension(:,:,:) :: p1Sp2_NB
    integer, save :: ibsize_print = OFF
#endif

#ifdef TRANSPOSE_WITHOUT_REARRANGEMENT
    integer :: ierror
#endif
    integer :: ibl1,ibl2,ibsize,ncache
    integer  :: id_sname = -1, id_sname2 = -1, id_sname3 = -1, id_sname4 = -1, id_sname1 = -1
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G ',id_sname)

    ncache = (m_CtrlP_cachesize()*1024)*3/4


    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(1) ',id_sname1)

#ifdef TRANSPOSE_WITHOUT_REARRANGEMENT
    if(ipri>=1) write(nfout,'("A CPP definition of " &
         & ,"TRANSPOSE_WITHOUT_REARRANGEMENT can not be set for MDDAVIDSON")')
    ierror = CPP_DEFINE_ERROR
#ifdef DEBUG_ERRORS
    call phase_error_wo_filename(ierror,nfout,line=__LINE__,modulefile=__FILE__)
#else
    call phase_error_wo_filename(ierror,nfout)
#endif
#endif

#ifdef MGS_DGEMM
    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
    if(ipri >= 2) then
       if(ibsize_print == OFF) then
          if(nblocksize_mgs_is_given) then
             write(nfout,'(" ! nblocksize_mgs_is_given")')
          else
             write(nfout,'(" ! nblocksize_mgs_is_given is false")')
          end if
          write(nfout,'( "! NB(=nblocksize_mgs) (mgs_phi2wf_each_k_G) = ",i8)') NB
          ibsize_print = ON
       end if
    end if
#endif

    call m_ESortho_mgs_alloc(ik)
    kimg_t_wk = 2
    if(k_symmetry(ik) == GAMMA) kimg_t_wk = 1
#ifdef MGS_DGEMM
    allocate(p1Sp2_NB(NB,NB,kimg_t_wk))
    allocate(wk_psi(mp_g1k(ik),NB,kimg) ) ; wk_psi = 0.0d0
    allocate(psi_t_dia(np_g1k_x,NB,kimg))
#endif

    dsize_psi = mp_g1k(ik)*NB*kimg
    if(mod_pot == VANDERBILT_TYPE) then
      if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
        i = 2
      else
        i = 4
      endif
      allocate( wk_bpri(i,max(mp_fs,nac_p),NB) )
      dsize_bpri = i*max(mp_fs,nac_p)*NB
    endif
    allocate( bpr_t_dia(np_fs,NB) )
    allocate( bpi_t_dia(np_fs,NB) )
    allocate( bpr_tw1_dia(nac_p,NB) )
    allocate( bpi_tw1_dia(nac_p,NB) )

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_r_3D(ista_k,iend_k,ik,fsr_l,bpr_t)    ! fsr_l -> bpr_t
          allocate(phifr_t(np_fs,np_e),phifi_t(1,1))
!!$          allocate(bsdr_lc(np_e,np_fs,ik:ik))
          call m_ES_F_transpose_r_3D(ik,    ik,    ik,bsdr_l,phifr_t) ! bsdr_l -> phifr_t
!!$          do i = 1, np_e
!!$             do j = ista_fs, iend_fs
!!$                bsdr_lc(i,j-ista_fs+1,ik) = bsdr_l(i,j,ik)
!!$             end do
!!$          end do
!!$          call m_ES_F_transpose_r_3D(ik,    ik,    ik,bsdr_lc,phifr_t) ! bsdr_l -> phifr_t
!!$#ifdef MGS_DGEMM
!!$          allocate(bpr_t_dia(np_fs_x,NB))
!!$#endif
       else
          call m_ES_F_transpose_r_3D(ista_k,iend_k,ik,fsr_l,bpr_t,fsi_l,bpi_t)  ! fs[ri]_l -> bp[ri]_t
          allocate(phifr_t(np_fs,np_e),phifi_t(np_fs,np_e))
!!$          allocate(bsdr_lc(np_e,np_fs,ik:ik),bsdi_lc(np_e,np_fs,ik:ik))
!!$          do i = 1, np_e
!!$             do j = ista_fs, iend_fs
!!$                bsdr_lc(i,j-ista_fs+1,ik) = bsdr_l(i,j,ik)
!!$                bsdi_lc(i,j-ista_fs+1,ik) = bsdi_l(i,j,ik)
!!$             end do
!!$          end do
          call m_ES_F_transpose_r_3D(ik,    ik,    ik,bsdr_l,phifr_t,bsdi_l,phifi_t)
!!$          call m_ES_F_transpose_r_3D(ik,    ik,    ik,bsdr_lc,phifr_t,bsdi_lc,phifi_t)
       end if
    end if

    allocate(phi_t(np_g1k_x,np_e,kimg))

    call m_ES_W_transpose_r_3D(ista_k,iend_k,ik,zaj_l,psi_t) ! zaj_l(ig,ie,ik,ri) -> psi_t(ig,ie,ri)
    call m_ES_W_transpose_r_3D(ik,    ik,    ik,phi_l,phi_t) ! phi_l(ig,ie,ik,ri) -> phi_t(ig,ie,ri)

    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname1)

    icount = (neg/NB+1)/nrank_e+1
!!$    allocate( p1Sp2_t2_NB(NB,NB,kimg_t_wk,icount) )
    allocate( p1Sp2_t1_NB(NB,NB,kimg_t_wk,icount) )

#ifdef MGS_DGEMM
    do i = 1, neg,NB
       NB_END = i + NB -1
       if( NB_END > neg ) NB_END = neg

       nbs = (i-1)/NB+1
       if( myrank_e == lrank(nbs)) then
          local_block = nbsn(nbs)
          L_NB_STA = nbsn_sta(local_block)
          L_NB_END = nbsn_end(local_block)
    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(2) ',id_sname2)
!diagonal
          do i1 = L_NB_STA, L_NB_END
             if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
                call WSW_t_g(ik,i1,mod_pot,fr,phi_t,np_g1k_x,np_e,kimg,phifr_t,phifi_t)
                if(dabs(fr-1.d0) > DELTA) &
                     & call normalize_bp_and_psi_t_g(ik,i1,fr,mod_pot &
                     &   , phi_t,np_g1k_x,np_e,kimg,phifr_t,phifi_t)
                !       |Phi(i)> = |Phi(i)> * fr,  <beta|Phi(i)> = <beta|Phi(i)> * fr
!x!!$          if(mod_pot == VANDERBILT_TYPE) call cp_bpr2bprtw(i1,L_NB_END)
             end if
             if(mod_pot == VANDERBILT_TYPE) call cp_bpr_bsdr2bprtw(i1,L_NB_END)
             !                  phif[ri]_t, bp[ri]_t -> bpr_tw1, bpr_tw2
             if(i1 == neg) cycle
             call cp_psi2psii_g(ik,i1) ! psi_t(:,i,:) -> psi_ir,psi_ii
             if(mod_pot == VANDERBILT_TYPE) &
                  & call cp_bpr2bpi_g(kimg_t_wk,i1,bpr_t,bpi_t) ! -> bp_ir, bp_ii
             call W1SW2_t_r_g(ik,i1,L_NB_END,mod_pot,phi_t,np_g1k_x,np_e,kimg) ! -> p1Sp2
             call modify_bp_and_psi_t_r_g(ik,i1,L_NB_END,mod_pot &
                  & ,phi_t,np_g1k_x,np_e,kimg,phifr_t,phifi_t)
             !                   psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
          end do   ! i1-loop
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
!lower

    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(3) ',id_sname3)
!!$    call cp_psi_bpri2dias_g(ik,i,kimg_t_wk,NB,psi_t_dia,np_g1k_x,NB,kimg &
!!$         &                 ,bpr_t_dia,bpi_t_dia,np_fs_x)
          call cp_psi_bpri2dias_g1_3D() ! -> wk_psi

       endif ! myrank_e == lrank(nbs)

       call make_diagonal()  ! wk_psi -> psi_t_dia, bp[ri]_t_dia, bp[ri]_tw1_dia

       icount = 0
       p1Sp2_t1_NB = 0.0d0

       do i2 = i+NB, neg,NB
          if(NB_END == neg) cycle
          NB_END2 = i2 + NB -1
          if( NB_END2 > neg ) NB_END2 = neg

          nbs = (i2-1)/NB+1
          if( myrank_e /= lrank(nbs)) cycle
          local_block = nbsn(nbs)
          L_NB_STA = nbsn_sta(local_block)
          L_NB_END = nbsn_end(local_block)
          if(mod_pot == VANDERBILT_TYPE) call cp_bpr_bsdr2bprtw(L_NB_STA,L_NB_END)
          !                                 phif[ri]_t, bp[ri]_t -> bpr_tw1, bpr_tw2

          call W1SW2_t_r_block_g(ik,i,L_NB_STA,p1Sp2_NB,NB_END,L_NB_END,kimg_t_wk &
               & , mod_pot,phi_t,psi_t_dia,np_g1k_x,np_e,kimg,bpr_tw1_dia,bpi_tw1_dia) ! -> p1Sp2
          !X!!$            & , mod_pot,phi_t,psi_t_dia,np_g1k_x,np_e,kimg,bpr_tw1(1,i),bpi_tw1(1,i)) ! -> p1Sp2
          icount = icount + 1

          do i_NB = 1, NB
             do iq = 1, kimg_t_wk
                do j_NB = 1,NB
                   p1Sp2_t1_NB(j_NB,i_NB,iq,icount) = p1Sp2_NB(j_NB,i_NB,iq)
          end do; enddo ; enddo
       end do ! i2-looop1

       if(NB_END /= neg) then
          if(nrank_g > 1 ) then
             ix =  NB*NB*kimg_t_wk*icount
             call mpi_allreduce(MPI_IN_PLACE,p1Sp2_t1_NB,ix &
                  & ,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          endif
       endif   !(NB_END /= neg )

       icount = 0
       do i2 = i+NB, neg,NB
          if(NB_END == neg) cycle
          NB_END2 = i2 + NB -1
          if( NB_END2 > neg ) NB_END2 = neg
          nbs = (i2-1)/NB+1
          if( myrank_e /= lrank(nbs)) cycle
          local_block = nbsn(nbs)
          L_NB_STA = nbsn_sta(local_block)
          L_NB_END = nbsn_end(local_block)
!modify2010
          icount =  icount+1
          do i_NB = 1, NB
             do iq = 1, kimg_t_wk
                do j_NB = 1, NB
                   p1Sp2_NB(j_NB,i_NB,iq) = p1Sp2_t1_NB(j_NB,i_NB,iq,icount)
          end do; enddo ; enddo

          call modify_bp_and_psi_t_r_blk_g(ik,i,L_NB_STA,p1Sp2_NB,NB_END,L_NB_END,kimg_t_wk &
               & , mod_pot,phi_t,psi_t_dia,np_g1k_x,np_e,kimg,phifr_t,phifi_t,bpr_t_dia,bpi_t_dia)
          ! phi_t, phifr_t, phifi_t, p1Sp2 -> phi_t, phifr_t, phifi_t
       end do   ! i2-loop
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname3)
    end do   ! i-loop
#else
    do i = 1, neg
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,i,mod_pot,fr,phi_t,np_g1k_x,np_e,kimg,phifr_t,phifi_t)
          if(dabs(fr-1.d0) > DELTA)  &
               & call normalize_bp_and_psi_t_g(ik,i,fr,mod_pot &
               &   , phi_t,np_g1k_x,np_e,kimg,phifr_t,phifi_t)
       end if
!!$       if(mod_pot == VANDERBILT_TYPE) call cp_bpr_bsdr2bprtw(i,neg)
!!$       if(i == neg) cycle
          call cp_psi2psii_g(ik,i) ! psi_t(:,i,:) -> psi_ir,psi_ii
       if(mod_pot == VANDERBILT_TYPE) &
            & call cp_bpr2bpi_g(kimg_t_wk,i,bpr_t,bpi_t) ! -> bp_ir, bp_ii
       call W1SW2_t_r_g(ik,i,neg,mod_pot,phi_t,np_g1k_x,np_e,kimg & ! ->p1Sp2
            & ,phifr_t,phifi_t)
       call modify_bp_and_psi_t_r_g(ik,i,neg,mod_pot &
            & ,phi_t,np_g1k_x,np_e,kimg,phifr_t,phifi_t)
       ! phi_t, psi_t, phifr_t,phifi_t,bpr_t, pbi_t, p1Sp2 -> phi_t, phifr_t,phifi_t

    end do
#endif

    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(4) ',id_sname4)
    call m_ES_W_transpose_back_r_3D(ik,ik,ik,phi_l,phi_t)  ! phi_t -> phi_l
    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_back_r_3D(ik,ik,ik,bsdr_l,phifr_t)
!!$          call m_ES_F_transpose_back_r_3D(ik,ik,ik,bsdr_lc,phifr_t)
!!$          bsdr_l(:,:,ik) = 0.d0
!!$          do i = 1, np_e
!!$             do j = ista_fs, iend_fs
!!$                bsdr_l(i,j,ik) = bsdr_lc(i,j-ista_fs+1,ik)
!!$             end do
!!$          end do
!!$          call mpi_allreduce(MPI_IN_PLACE,bsdr_l(1,1,ik),np_e*nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
!!$          deallocate(bsdr_lc)
       else
          call m_ES_F_transpose_back_r_3D(ik,ik,ik,bsdr_l,phifr_t,bsdi_l,phifi_t)
!!$          call m_ES_F_transpose_back_r_3D(ik,ik,ik,bsdr_lc,phifr_t,bsdi_lc,phifi_t)
!!$          bsdr_l(:,:,ik) = 0.d0
!!$          bsdi_l(:,:,ik) = 0.d0
!!$          do i = 1, np_e
!!$             do j = ista_fs, iend_fs
!!$                bsdr_l(i,j,ik) = bsdr_lc(i,j-ista_fs+1,ik)
!!$                bsdi_l(i,j,ik) = bsdi_lc(i,j-ista_fs+1,ik)
!!$             end do
!!$          end do
!!$          call mpi_allreduce(MPI_IN_PLACE,bsdr_l(1,1,ik),np_e*nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
!!$          call mpi_allreduce(MPI_IN_PLACE,bsdi_l(1,1,ik),np_e*nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
!!$          deallocate(bsdr_lc,bsdi_lc)
       end if
    end if

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          deallocate(phifr_t)
       else
          deallocate(phifr_t,phifi_t)
       end if
    end if

    deallocate(p1Sp2_t1_NB)
    deallocate(phi_t)
    deallocate(bpi_tw1_dia,bpr_tw1_dia)
    deallocate(bpi_t_dia,bpr_t_dia)
    if(mod_pot == VANDERBILT_TYPE) deallocate(wk_bpri)
#ifdef MGS_DGEMM
    deallocate(psi_t_dia)
    deallocate(wk_psi)
    deallocate(p1Sp2_NB)
#endif

    call m_ESortho_mgs_dealloc
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname4)
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  contains
!!$#ifdef MGS_DGEMM
!!$    subroutine cp_bpr2bprtw(i1,i2)
!!$      integer, intent(in) :: i1,i2
!!$      integer :: i, ia, p, q, i_last
!!$      i_last = i2
!!$      if(k_symmetry(ik) == GAMMA) then
!!$         do i = i1, i_last
!!$            do ia = 1, nac_p
!!$               p = nlmta1_p(ia);     q = nlmta2_p(ia)
!!$               bpr_tw1(ia,i) = bpr_t(p,i)
!!$               bpr_tw2(ia,i) = bpr_t(q,i)
!!$            end do
!!$         end do
!!$      else
!!$         do i = i1, i_last
!!$            do ia = 1, nac_p
!!$               p = nlmta1_p(ia);     q = nlmta2_p(ia)
!!$               bpr_tw1(ia,i) = bpr_t(p,i)
!!$               bpr_tw2(ia,i) = bpr_t(q,i)
!!$               bpi_tw1(ia,i) = bpi_t(p,i)
!!$               bpi_tw2(ia,i) = bpi_t(q,i)
!!$            end do
!!$         end do
!!$      end if
!!$    end subroutine cp_bpr2bprtw
!!$#endif

#ifdef MGS_DGEMM
    subroutine make_diagonal()
      integer :: ri,ix,iy

      if(nrank_e > 1) then
         call mpi_bcast(wk_psi,dsize_psi,mpi_double_precision,lrank(nbs),mpi_kg_world,ierr)
      end if
      do ri = 1, kimg
         do iy = 1, NB
            do ix = 1, np_g1k(ik)
               psi_t_dia(ix,iy,ri) = wk_psi(ix,iy,ri)
            end do
         end do
      end do

      if(mod_pot == VANDERBILT_TYPE) then
         if(nrank_e > 1) then
            call mpi_bcast(wk_bpri,dsize_bpri &
                 & ,mpi_double_precision,lrank(nbs),mpi_kg_world,ierr)
         end if
         if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
            do iy = 1, NB
               do ix = 1, np_fs
                  bpr_t_dia(ix,iy) = wk_bpri(1,ix,iy)
               enddo
               do ix = 1, nac_p
                  bpr_tw1_dia(ix,iy) = wk_bpri(2,ix,iy)
               enddo
            enddo
         else
            do iy = 1, NB
               do ix = 1, np_fs
                  bpr_t_dia(ix,iy) = wk_bpri(1,ix,iy)
                  bpi_t_dia(ix,iy) = wk_bpri(2,ix,iy)
               enddo
               do ix = 1, nac_p
                  bpr_tw1_dia(ix,iy) = wk_bpri(3,ix,iy)
                  bpi_tw1_dia(ix,iy) = wk_bpri(4,ix,iy)
               enddo
            enddo
         end if
      end if
    end subroutine make_diagonal

    subroutine cp_psi_bpri2dias_g1_3D()
      integer :: ri,i1,ix,iy

      do ri = 1, kimg
         do i1 = L_NB_STA, L_NB_END
            iy = i1-L_NB_STA+1
            do ix = 1, np_g1k(ik)
               wk_psi(ix,iy,ri) = psi_t(ix,i1,ri)
            enddo
         enddo
      enddo

      if(mod_pot == VANDERBILT_TYPE) then
         if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
            do i1 = L_NB_STA, L_NB_END
               iy = i1-L_NB_STA+1
               do ix = 1, np_fs
                  wk_bpri(1,ix,iy) = bpr_t(ix,i1)
               enddo
               do ix = 1, nac_p
                  wk_bpri(2,ix,iy) = bpr_tw1(ix,i1)
               enddo
            enddo
         else
            do i1 = L_NB_STA, L_NB_END
               iy = i1-L_NB_STA+1
               do ix = 1, np_fs
                  wk_bpri(1,ix,iy) = bpr_t(ix,i1)
                  wk_bpri(2,ix,iy) = bpi_t(ix,i1)
               enddo
               do ix = 1, nac_p
                  wk_bpri(3,ix,iy) = bpr_tw1(ix,i1)
                  wk_bpri(4,ix,iy) = bpi_tw1(ix,i1)
               enddo
            enddo
         end if
      end if
    end subroutine cp_psi_bpri2dias_g1_3D

    subroutine cp_bpr_bsdr2bprtw(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: i, ia, p, q, i_last
#ifdef MGS_DGEMM_DEBUG
      i_last = i1
#else
      i_last = i2
#endif
      if(k_symmetry(ik) == GAMMA) then
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = phifr_t(q,i)
            end do
         end do
      else
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = phifr_t(q,i)
               bpi_tw1(ia,i) = bpi_t(p,i)
               bpi_tw2(ia,i) = phifi_t(q,i)
!!$               bpr_tw1(ia,i) = phifr_t(p,i)
!!$               bpr_tw2(ia,i) = bpr_t(q,i)
!!$               bpi_tw1(ia,i) = phifi_t(p,i)
!!$               bpi_tw2(ia,i) = bpi_t(q,i)
            end do
         end do
      end if
    end subroutine cp_bpr_bsdr2bprtw
#endif
  end subroutine mgs_phi2wf_each_k_G

!******************** modified by RIST_16
  subroutine WSW_for_overlap_t_g(ik,ito,jto,fr,fi,bpr_t,bpi_t)
    use m_Electronic_Structure, only : fsr_l, fsi_l
    integer, intent(in)        :: ik,ito,jto
    real(kind=DP), intent(out) :: fr,fi
    real(kind=DP),intent(in),optional :: bpr_t(:,:),bpi_t(:,:)
    real(kind=DP),allocatable :: bpr(:,:),bpi(:,:)

    integer              :: ia,p,q, i, ig1,ib
    integer :: id_sname = -1
    allocate(bpr(np_fs,neg))
    if(.not.(k_symmetry(ik)==GAMMA.and.kimg == 2)) allocate(bpi(np_fs,neg))
    bpr=0.d0
    do ib = 1,np_e
       do i = 1, np_fs
          !bpr  (i,ib) = fsr_l(ib-ista_e+1,i,ik)
          bpr  (i,neg_g(ib)) = bpr_t(ib,i)
       end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,bpr,  np_fs*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
    if(allocated(bpi)) then
      bpi=0.d0
      do ib = 1,np_e
         do i = 1, np_fs
            !bpi  (i,ib) = fsi_l(ib-ista_e+1,i,ik)
            bpi  (i,neg_g(ib)) = bpi_t(ib,i)
         end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,bpi,  np_fs*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
    endif

    fr = 0.d0
    fi = 0.d0
    if(k_symmetry(ik) == GAMMA .and. kimg == 2) then
       do ia = 1, nac_p
          p = nlmta1_p(ia);     q = nlmta2_p(ia)
          fr = fr+fqwei_p(ia)*(bpr(p,jto)*bpr(q,ito))
       end do
    else
       do ia = 1, nac_p
          p = nlmta1_p(ia);     q = nlmta2_p(ia)
          fr = fr+fqwei_p(ia)*(bpr(p,jto)*bpr(q,ito)+bpi(p,jto)*bpi(q,ito))
          fi = fi+fqwei_p(ia)*(bpr(p,jto)*bpi(q,ito)-bpi(p,jto)*bpr(q,ito))
       end do
    end if
    !call mpi_allreduce(MPI_IN_PLACE,fr,1,mpi_double_precision,MPI_SUM,mpi_ke_world,ierr)
    !if(allocated(bpi)) &
    !&  call mpi_allreduce(MPI_IN_PLACE,fi,1,mpi_double_precision,MPI_SUM,mpi_ke_world,ierr)
    deallocate(bpr)
    if(allocated(bpi)) deallocate(bpi)
  end subroutine WSW_for_overlap_t_g
end module m_ES_ortho
