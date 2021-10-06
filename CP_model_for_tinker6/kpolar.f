c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c
      subroutine kpolar
! qtw >> in
      use chgpen, only : ecp
! qtw << end
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kpolr.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'usolve.i'
      integer i,j,k
      integer npg,next
      integer pg(maxval)
      real*8 pol,thl
      real*8 sixth
      logical header
      character*20 keyword
      character*120 record
      character*120 string
! qtw >> in
      integer :: n_valence
      real*8  :: a_jp
      real*8  :: b_jp
      real*8  :: s_lo
      real*8  :: s_hi
      logical :: chgpen_header
! qtw << out
c
c
c     process keywords containing polarizability parameters
c
      header = .true.
      chgpen_header = .true.  ! qtw 
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:120)
            read (string,*,err=10,end=10)  pol,thl,(pg(j),j=1,maxval)
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,5x,'Atom Type',11x,'Alpha',8x,
     &                       'Damp',5x,'Group Atom Types'/)
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  athl(k) = thl
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 30
                     end if
                  end do
   30             continue
                  if (.not. silent) then
                     write (iout,40)  k,pol,thl,(pg(j),j=1,npg)
   40                format (4x,i6,10x,f10.3,2x,f10.3,7x,20i5)
                  end if
               else
                  write (iout,50)
   50             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
! qtw >> in
c
c     charge penetration parameters
c
         else if (keyword(1:10) .eq. 'CHGPENPRM ') then
            k = 0
            n_valence = 0
            a_jp = 0.0D0
            b_jp = 0.0D0

            call getnumb (record,k,next)
            string = record(next:120)
            ! read type N_val_JP alpha_JP beta_JP
            read (string,*,err=200,end=200) n_valence,a_jp, b_jp
  200       continue
            if (k .gt. 0) then
               if (chgpen_header) then
                  chgpen_header = .false.
                  write (iout,210)
  210             format (/,' Additional Charge',
     &                        ' Penetration Parameters:',
     &                    /,6x,'Atom',15x,'N_val',7x,'Alpha',7x,'Beta'/)
               end if

               if (k .le. maxtyp) then
                  ecp%n_val(k)  = n_valence
                  ecp%alp_jp(k) = a_jp
                  ecp%bet_jp(k) = b_jp

                  write (iout,220) k,n_valence,a_jp,b_jp
  220             format (4x,i6,10x,i10,2x,f10.3,2x,f10.3)
               else
                  write (iout,230)
  230             format (/,' KPOLAR  --  Too many',
     &                       ' Charge Penetration Parameters')
                  abort = .true.
               end if

            end if

         else if (keyword(1:14) .eq. 'CHGPEN-SWITCH ') then
            s_lo = 0.0D0
            s_hi = 0.0D0

            string = record(next:120)
            ! read switch_lo and switch_hi
            read (string,*,err=446,end=446) s_lo, s_hi
  446       continue
            ecp%switch_lo = s_lo
            ecp%switch_hi = s_hi

! qtw >> out
         end if
      end do
c
c     find and store the atomic dipole polarizability parameters
c
      do i = 1, n
         polarity(i) = polr(type(i))
         thole(i) = athl(type(i))
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      chgpen_header = .true.  ! qtw
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
               read (string,*,err=60,end=60)  pol,thl
   60          continue
               if (header) then
                  header = .false.
                  write (iout,70)
   70             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',8x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,80)  k,pol,thl
   80             format (4x,i6,10x,f10.3,2x,f10.3)
               end if
               polarity(k) = pol
               thole(k) = thl
            end if
! qtw >> in
c
c     charge penetration parameters
c
         else if (keyword(1:10) .eq. 'CHGPENPRM ') then
            k = 0
            n_valence = 0
            a_jp = 0.0D0
            b_jp = 0.0D0

            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:120)
            ! read type N_val_JP alpha_JP beta_JP
               read (string,*,err=240,end=240) n_valence,a_jp,b_jp
  240          continue
               if (chgpen_header) then
                  chgpen_header = .false.
                  write (iout,250)
  250             format (/,' Additional Charge Penetration Parameters',
     &                       ' for Specific Atoms :',
     &                    /,6x,'Atom',15x,'N_val',7x,'Alpha',7x,'Beta'/)
               end if
               write (iout,260) k,n_valence,a_jp,b_jp
  260          format (4x,i6,10x,i10,2x,f10.3,2x,f10.3)

               ecp%n_val(k)  = n_valence
               ecp%alp_jp(k) = a_jp
               ecp%bet_jp(k) = b_jp
            end if

         else if (keyword(1:14) .eq. 'CHGPEN-SWITCH ') then
            s_lo = 0.0D0
            s_hi = 0.0D0

            string = record(next:120)
            ! read switch_lo and switch_hi
            read (string,*,err=447,end=447) s_lo, s_hi
  447       continue
            ecp%switch_lo = s_lo
            ecp%switch_hi = s_hi

! qtw >> out
         end if
      end do
c
c     remove zero and undefined polarizable sites from the list
c
      npole = 0
      npolar = 0
      do i = 1, n
         if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
            npole = npole + 1
            ipole(npole) = i
            pollist(i) = npole
            zaxis(npole) = zaxis(i)
            xaxis(npole) = xaxis(i)
            yaxis(npole) = yaxis(i)
            polaxe(npole) = polaxe(i)
            do k = 1, maxpole
               pole(k,npole) = pole(k,i)
            end do
            if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
            polarity(npole) = polarity(i)
            thole(npole) = thole(i)
         end if
      end do
c
c     set the values used in the scaling of the polarizability
c
      sixth = 1.0d0 / 6.0d0
      do i = 1, npole
         if (thole(i) .eq. 0.0d0) then
            pdamp(i) = 0.0d0
         else
            pdamp(i) = polarity(i)**sixth
         end if
      end do
c
c     assign polarization group connectivity of each atom
c
      call polargrp
c
c     test multipoles at chiral sites and invert if necessary
c
      call chkpole
c
c     turn off polarizable multipole potential if it is not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (npolar .eq. 0)  use_polar = .false.
c
c     perform dynamic allocation of some pointer arrays
c
      if (use_polar) then
         if (associated(mindex))  deallocate (mindex)
         if (associated(minv))  deallocate (minv)
         allocate (mindex(npole))
         allocate (minv(3*maxulst*npole))
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kpolr.i'
      include 'mpole.i'
      include 'polgrp.i'
      integer maxlist,maxkeep
      parameter (maxkeep=100)
      parameter (maxlist=1000)
      integer i,j,k
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nlist,nkeep
      integer keep(maxkeep)
      integer list(maxlist)
      integer, allocatable :: mask(:)
      logical done
c
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         it = type(i)
         do j = 1, n12(i)
            jj = i12(j,i)
            jt = type(jj)
            do k = 1, maxval
               kk = pgrp(k,it)
               if (kk .eq. 0)  goto 20
               if (pgrp(k,it) .eq. jt) then
                  np11(i) = np11(i) + 1
                  if (np11(i) .le. maxp11) then
                     ip11(np11(i),i) = jj
                  else
                     write (iout,10)
   10                format (/,' POLARGRP  --  Too many Atoms',
     &                          ' in Polarization Group')
                     abort = .true.
                     goto 30
                  end if
               end if
            end do
   20       continue
         end do
      end do
   30 continue
c
c     perform dynamic allocation of some local arrays
c
      allocate (mask(n))
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               mask(jj) = i
            end if
         end do
         do while (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (mask(kk) .ne. i) then
                     np11(i) = np11(i) + 1
                     if (np11(i) .le. maxp11) then
                        ip11(np11(i),i) = kk
                     else
                        write (iout,40)
   40                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 50
                     end if
                     mask(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
   50 continue
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp12) then
            np12(i) = nlist
            do j = 1, nlist
               ip12(j,i) = list(j)
            end do
         else
            write (iout,60)
   60       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-2 Polarization Group')
            abort = .true.
            goto 70
         end if
      end do
   70 continue
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp13) then
            np13(i) = nlist
            do j = 1, nlist
               ip13(j,i) = list(j)
            end do
         else
            write (iout,80)
   80       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-3 Polarization Group')
            abort = .true.
            goto 90
         end if
      end do
   90 continue
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp14) then
            np14(i) = nlist
            do j = 1, nlist
               ip14(j,i) = list(j)
            end do
         else
            write (iout,100)
  100       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-4 Polarization Group')
            abort = .true.
            goto 110
         end if
      end do
  110 continue
c
c     perform deallocation of some local arrays
c
      deallocate (mask)
      return
      end
