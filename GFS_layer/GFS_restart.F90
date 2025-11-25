module physics_restart_layer

  use machine,                   only: kind_phys
  use IPD_typedefs,              only: IPD_restart_type
  use physics_abstraction_layer, only: control_type,  statein_type,  &
                                       stateout_type, sfcprop_type,  &
                                       coupling_type, grid_type,     &
                                       statemid_type, cldprop_type,  &
                                       radtend_type,  intdiag_type,  &
                                       init_type   

  public restart_populate

  CONTAINS
!*******************************************************************************************

!---------------------
! GFS_restart_populate
!---------------------
  subroutine restart_populate (IPD_Restart, Model, Statein, Stateout, Sfcprop, &
                               Coupling, Grid, Statemid, Cldprop, Radtend, Diag, Init_parm)
!----------------------------------------------------------------------------------------!
!   IPD_METADATA                                                                         !
!     IPD_Restart%num2d          [int*4  ]  number of 2D variables to output             !
!     IPD_Restart%num3d          [int*4  ]  number of 3D variables to output             !
!     IPD_Restart%name2d         [char=32]  variable name in restart file                !
!     IPD_Restart%name3d         [char=32]  variable name in restart file                !
!     IPD_Restart%fld2d(:,:,:)   [real*8 ]  pointer to 2D data (im,nblks,MAX_RSTRT)      !
!     IPD_Restart%fld3d(:,:,:,:) [real*8 ]  pointer to 3D data (im,levs,nblks,MAX_RSTRT) !
!----------------------------------------------------------------------------------------!
    type(IPD_restart_type),  intent(inout) :: IPD_Restart
    type(control_type),         intent(in)    :: Model
    type(statein_type),         intent(in)    :: Statein(:)
    type(stateout_type),        intent(in)    :: Stateout(:)
    type(sfcprop_type),         intent(in)    :: Sfcprop(:)
    type(coupling_type),        intent(in)    :: Coupling(:)
    type(grid_type),            intent(in)    :: Grid(:)
    type(statemid_type),        intent(in)    :: Statemid(:)
    type(cldprop_type),         intent(in)    :: Cldprop(:)
    type(radtend_type),         intent(in)    :: Radtend(:)
    type(intdiag_type),         intent(in)    :: Diag(:)
    type(init_type),            intent(in)    :: Init_parm

    !--- local variables
    integer :: nblks, num, nb, max_rstrt, offset
    character(len=2) :: c2 = ''
    
    nblks = size(Init_parm%blksz)
    max_rstrt = size(IPD_Restart%name2d)

    !TODO lmh 14 jan 2020
    ! The MLO variables should really be saved in sfc_restart
    !  and not phy_restart.
    IPD_Restart%num2d = 5 + 10 + Model%ntot2d + Model%nctp
    IPD_Restart%num3d = Model%ntot3d

    allocate (IPD_Restart%name2d(IPD_Restart%num2d))
    allocate (IPD_Restart%name3d(IPD_Restart%num3d))
    allocate (IPD_Restart%data(nblks,max(IPD_Restart%num2d,IPD_Restart%num3d)))

    IPD_Restart%name2d(:) = ' '
    IPD_Restart%name3d(:) = ' '

    !--- Cldprop variables
    IPD_Restart%name2d(1) = 'cv'
    IPD_Restart%name2d(2) = 'cvt'
    IPD_Restart%name2d(3) = 'cvb'
    do nb = 1,nblks
      IPD_Restart%data(nb,1)%var2p => Cldprop(nb)%cv(:)
      IPD_Restart%data(nb,2)%var2p => Cldprop(nb)%cvt(:)
      IPD_Restart%data(nb,3)%var2p => Cldprop(nb)%cvb(:)
    enddo

    !--- Mixed-layer ocean variables
    offset = 3
    IPD_Restart%name2d(1+offset) = 'ts_som'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%ts_som(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'tsclim'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%tsclim(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'mldclim'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%mldclim(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'ts_clim_iano'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%ts_clim_iano(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'tml'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%tml(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'tml0'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%tml0(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'mld'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%mld(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'mld0'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%mld0(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'huml'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%huml(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'hvml'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%hvml(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'tmoml'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%tmoml(:)
    enddo

    offset = offset + 1
    IPD_Restart%name2d(1+offset) = 'tmoml0'
    do nb = 1,nblks
      IPD_Restart%data(nb,1+offset)%var2p => Sfcprop(nb)%tmoml0(:)
    enddo

    !TODO lmh 14 jan 2020
    ! Most of the phy_restart variables are redundant with the
    !  tracers saved in the dynamics, and are not needed.
    !--- phy_f2d variables
    offset = offset + 1
    do num = 1,Model%ntot2d
       !--- set the variable name
      write(c2,'(i2.2)') num
      IPD_Restart%name2d(num+offset) = 'phy_f2d_'//c2
      do nb = 1,nblks
        IPD_Restart%data(nb,num+offset)%var2p => Statemid(nb)%phy_f2d(:,num)
      enddo
    enddo

    !--- phy_fctd variables
    offset = offset + Model%ntot2d
    do num = 1, Model%nctp
       !--- set the variable name
      write(c2,'(i2.2)') num
      IPD_Restart%name2d(num+offset) = 'phy_fctd_'//c2
      do nb = 1,nblks
        IPD_Restart%data(nb,num+offset)%var2p => Statemid(nb)%phy_fctd(:,num)
      enddo
    enddo

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
       !--- set the variable name
      write(c2,'(i2.2)') num
      IPD_Restart%name3d(num) = 'phy_f3d_'//c2
      do nb = 1,nblks
        IPD_Restart%data(nb,num)%var3p => Statemid(nb)%phy_f3d(:,:,num)
      enddo
    enddo

  end subroutine restart_populate

end module physics_restart_layer
