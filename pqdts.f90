!pqdts: parallel qauntum detector tomography solver
!MIT License
!
!Copyright (c) 2024 Paderborn University, Paderborn Center for Parallel Computing, Robert Schade
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.

module basics
  use omp_lib, only: omp_get_wtime
  use iso_fortran_env, only: int8, int32, int64, real64
#ifdef MPI
  use mpi_f08, ONLY: MPI_request, MPI_Barrier, MPI_COMM_WORLD, MPI_allreduce, MPI_DOUBLE, MPI_SUM, &
                     MPI_Status, MPI_STATUSES_IGNORE, MPI_RECV, MPI_MIN, mpi_in_place, mpi_long, mpi_max, MPI_Finalize
#endif

  implicit none

  integer(kind=int64) :: m, n, d, MB, MB1, MB2, NP1, NP2
  integer(kind=int64) :: Fl, Pl

  type :: dpline
    integer(kind=int32) :: c = 0
    logical :: isalloc = .False.
    real(kind=real64), pointer, dimension(:) :: f
  end type
  type :: intline
    integer(kind=int32) :: c = 0
    integer(kind=int32), pointer, dimension(:) :: f
  end type

#ifdef MPI
  !stuff for MPI
  integer(kind=int32) :: ierror
  integer(kind=int32) :: ierr
  TYPE(MPI_Request), pointer, dimension(:) :: requests
  TYPE(dpline), pointer, dimension(:) :: requestsbuf
  logical :: tbarrier
#endif
  integer(kind=int32) :: rank
  integer(kind=int32) :: num_proc
  integer(kind=int64) :: localrows
  integer(kind=int64) :: localrows0 !start of block of rows
  integer(kind=int64) :: localrows1
  type(intline), dimension(:), pointer :: block_red
  type(intline), dimension(:), pointer :: block_red_steps
  integer(int32) :: blocks_red

  !timing
  real(kind=real64), pointer, dimension(:) :: timingsum
  real(kind=real64), pointer, dimension(:) :: timingt
  integer(kind=int32), pointer, dimension(:) :: timingc
  integer(kind=int32) :: timingid
  character(len=20), pointer, dimension(:) :: timingn
  logical :: time_io = .false.

  !traffic
  real(kind=real64), pointer, dimension(:) :: dsend
  real(kind=real64), pointer, dimension(:) :: drecv

  !F
  integer(kind=int64), pointer, dimension(:) :: Fi, Fj
  real(kind=real64), pointer, dimension(:) :: Fd
  integer(kind=int64), pointer, dimension(:) :: Frow0
  integer(kind=int64), pointer, dimension(:) :: Frow1
  integer(kind=int64), pointer, dimension(:) :: Flocalrow0
  integer(kind=int64), pointer, dimension(:) :: Flocalrow1
  type(dpline), dimension(:), pointer :: Frows
  !P
  integer(kind=int64), pointer, dimension(:) :: Pi, Pj
  real(kind=real64), pointer, dimension(:) :: Pd

  !auxiliary variables
  real(kind=real64), pointer, dimension(:, :) :: P
  real(kind=real64), pointer, dimension(:, :) :: O
  real(kind=real64), pointer, dimension(:, :) :: Xtmp, Xp, sn !,btmp,btmp2,btmp3
  integer(kind=int8), pointer, dimension(:, :) :: ik
  integer(kind=int32), pointer, dimension(:) :: maxr

  !mode setting
  integer(kind=int32) :: mode
  logical :: tfullspace

  !constants
  real(kind=real64) :: M_PI = 4.D0*DATAN(1.D0)
  real(kind=real64) :: gam = 0.0D0
contains
  subroutine tstart(n)
    implicit none
    character(*), intent(in) :: n
    logical :: found
    integer(kind=int32) :: i, it
    real(kind=real64) :: t

    t = omp_get_wtime()
    found = .false.
    do i = 1, 50
      if (trim(n) .eq. timingn(i)) then
        it = i
        found = .true.
        exit
      end if
    end do

    if (.not. found) then
      timingid = timingid + 1
      it = timingid
    end if

    if (timingt(it) .ne. 0) then
      print *, "timer ", it, "is already running"
      stop
    end if
    timingt(it) = t
    timingc(it) = timingc(it) + 1
    timingn(it) = trim(n)
  end subroutine

  subroutine tstop(n)
    implicit none
    character(*), intent(in) :: n
    logical :: found
    integer(kind=int32) :: i, it
    real(kind=real64) :: t
    t = omp_get_wtime()

    found = .false.
    do i = 1, 50
      if (trim(n) .eq. timingn(i)) then
        it = i
        found = .true.
        exit
      end if
    end do

    if (.not. found) then
      print *, "timer ", trim(n), "has not been started"
      stop
    end if
    timingsum(it) = timingsum(it) + t - timingt(it)
    timingt(it) = 0
  end subroutine

  subroutine tprint(iter)
    implicit none
    integer(kind=int32), intent(in) :: iter
    integer(kind=int32) :: i
    Character(len=255) :: str
    logical :: file_exists

    if (.not. time_io) then
      write (str, *) "timing_", rank, ".out"
      str = adjustl(str)
      INQUIRE (FILE=trim(str), EXIST=file_exists)

      if (file_exists) then
        if (iter .eq. 0) then
          open (43, file=trim(str), status="replace", action="write")
        else
          open (43, file=trim(str), status="old", position="append", action="write")
        end if
      else
        open (43, file=trim(str), status='new', action="write")
      end if
      time_io = .true.
    end if
    do i = 1, timingid
      write (43, *) iter, timingn(i), timingc(i), timingsum(i)
    end do
    call flush (43)

  end subroutine

  subroutine fprob(x, mean, prob)
    implicit none
    real(kind=real64), intent(in) :: x
    real(kind=real64), intent(in) :: mean
    real(kind=real64), intent(out) :: prob
    real(kind=real64) :: y, lf

    if (mean .lt. 50) then
      if (mean .eq. 0 .and. x .eq. 0) then
        prob = 1.0D0
      else
        lf = LGAMMA(x + 1.0D0)
        prob = exp(log(mean)*x - lf - mean)
      end if
    else
      y = (x - mean)/sqrt(mean)
      prob = exp(-0.5D0*y**2)/sqrt(mean)/sqrt(2*M_PI)
    end if
  end subroutine

  subroutine Bdiag(BD)
    !diagonal of Hessian H_ii
    implicit none
    real(kind=real64), intent(inout), dimension(:, :), pointer :: BD
    integer(kind=int64) :: i, j
    call tstart("Bdiag")
    call pset(m, n, 0.0D0, BD)

    !$OMP parallel DO private(i) schedule(static,1)
    do j = 1, N
      do i = 1, D
        if (Flocalrow0(i) .ne. -1) then
          BD(Flocalrow0(i):Flocalrow1(i), j) = BD(Flocalrow0(i):Flocalrow1(i), j) + 2.0*Frows(i)%f**2
        end if
      end do
    end do
    !$OMP end parallel DO
    if (gam .ne. 0) then
      call tstart("Biag_gam")
      !$OMP parallel DO schedule(static,1)
      do i = 1, N
        BD(:localrows - 1, i) = BD(:localrows - 1, i) + 2.0D0*gam
        BD(2:, i) = BD(2:, i) + 2.0D0*gam
      end do
      !$OMP end parallel DO
      call tstop("Biag_gam")
    end if
    if (.not. tfullspace) then
      call pscale(m, n, 2.0D0, BD)
!      !$OMP parallel DO schedule(static,1)
!      do j = 1, localrows
!        BD(j, maxr(j)) = 0
!      end do
!      !$OMP end parallel DO
    end if
    call tstop("Bdiag")
  end subroutine

  subroutine BF(dd, BX)
    !Hessian times vector: BX=H(X)*dd
    implicit none
    real(kind=real64), intent(inout), dimension(:, :), pointer :: dd, BX
    integer(kind=int64) :: i, j

    !call pset(m,n,0.0D0,btmp)
    !call dl(btmp,btmp2)
    !call dl(dd,btmp3)
    !BX=btmp3-btmp2
    !return

    call tstart("BF")
    call tstart("BF_prep")
    if (.not. tfullspace) then
      call pcopy(m, n, dd, bx)
      !$OMP parallel DO schedule(static,1)
      do j = 1, localrows
        Bx(j, maxr(j)) = -(sum(dd(j, :)) - dd(j, maxr(j)))
      end do
      !$OMP end parallel DO
    end if
    call tstop("BF_prep")

    call tstart("BF_O")
    if (.not. tfullspace) then
      call redO(Bx, .true., O)
    else
      call redO(dd, .true., O)
    end if
    call pset(m, n, 0.0D0, BX)
    call tstop("BF_O")

    call tstart("BF_BX")
    !$OMP parallel DO private(i) schedule(static,1)
    do j = 1, N
      do i = 1, D
        if (Flocalrow0(i) .ne. -1) then
          BX(Flocalrow0(i):Flocalrow1(i), j) = BX(Flocalrow0(i):Flocalrow1(i), j) + 2.0*Frows(i)%f*(O(j, i))
        end if
      end do
    end do
    !$OMP end parallel DO
    call tstop("BF_BX")

    if (gam .ne. 0) then
      call tstart("BF_gam")
      if (.not. tfullspace) then
        !$OMP parallel DO schedule(static,1) private(j)
        do i = 1, N
          do j = 1, localrows
            if (j .lt. localrows) then
              BX(j, i) = BX(j, i) + 2.0D0*gam*(dd(j, i) - dd(j + 1, i))
              if (i .eq. maxr(j)) then
                BX(j, i) = BX(j, i) + 2.0D0*gam*(-dd(j, i) - (sum(dd(j, :)) - dd(j, maxr(j))))
              end if
              if (i .eq. maxr(j + 1)) then
                BX(j, i) = BX(j, i) + 2.0D0*gam*(+dd(j + 1, i) + (sum(dd(j + 1, :)) - dd(j + 1, maxr(j + 1))))
              end if
            end if
            if (j .gt. 1) then
              BX(j, i) = BX(j, i) - 2.0D0*gam*(dd(j - 1, i) - dd(j, i))
              if (i .eq. maxr(j - 1)) then
                BX(j, i) = BX(j, i) - 2.0D0*gam*(-dd(j - 1, i) - (sum(dd(j - 1, :)) - dd(j - 1, maxr(j - 1))))
              end if
              if (i .eq. maxr(j)) then
                BX(j, i) = BX(j, i) - 2.0D0*gam*(+dd(j, i) + (sum(dd(j, :)) - dd(j, maxr(j))))
              end if
            end if
          end do
        end do
        !$OMP end parallel DO
      else
        !$OMP parallel DO schedule(static,1) private(j)
        do i = 1, N
          do j = 1, localrows
            if (j .lt. localrows) then
              BX(j, i) = BX(j, i) + 2.0D0*gam*(dd(j, i) - dd(j + 1, i))
            end if
            if (j .gt. 1) then
              BX(j, i) = BX(j, i) - 2.0D0*gam*(dd(j - 1, i) - dd(j, i))
            end if
          end do
        end do
        !$OMP end parallel DO
      end if
      call tstop("BF_gam")
    end if
    if (.not. tfullspace) then
      !$OMP parallel DO schedule(static,1)
      do j = 1, localrows
        BX(j, :) = BX(j, :) - BX(j, maxr(j))
      end do
      !$OMP end parallel DO
    end if
    call tstop("BF")
  end subroutine

  subroutine Lalpha_min(X, S, tproj, alpha, dlalpha)
    implicit none
    real(kind=real64), intent(inout), dimension(:, :) :: X, S
    logical, intent(in) :: tproj
    real(kind=real64), intent(out) :: alpha, dlalpha
    real(kind=real64) :: df, dx, Lm, Lp, L0, alphamin
    integer(kind=int64) :: i

    call tstart("Lalpha_min")
    !backtracking
    dx = 1e-5
    alpha = 0.0D0
    call pgaxpy2(m, n, 1.0D0, x, alpha, s, xtmp)
    if (tproj) call doproj(Xtmp)
    call L(Xtmp, L0)

    call pgaxpy2(m, n, 1.0D0, x, alpha + dx, s, xtmp)
    if (tproj) call doproj(Xtmp)
    call L(Xtmp, Lp)
    df = (Lp - L0)/(dx)
    dlalpha = df
    if (mode .eq. 2) then
      alpha = 1.0D0
    else
      alpha = 100.0D0
    end if

    lm = huge(lm)
    alphamin = 0

    do i = 1, 100
      call pgaxpy2(m, n, 1.0D0, x, alpha, s, xtmp)
      if (tproj) call doproj(Xtmp)
      call L(Xtmp, Lp)
      if (lp .lt. lm) then
        lm = lp
        alphamin = alpha
      end if
      if (rank .eq. 0) print *, "line search", i, "L=", Lp, "alpha=", alpha, "df(0)=", df
      if (Lp .gt. L0 + 0.1D0*alpha*df) then
        alpha = alpha*0.75D0
      else
        exit
      end if
    end do
    alpha = alphamin
    if (lm .gt. l0) alpha = 0
    call tstop("Lalpha_min")
  end subroutine

  subroutine Lalpha_min2(X, S, tproj, alpha, dlalpha)
    implicit none
    real(kind=real64), intent(inout), dimension(:, :) :: X, S
    logical, intent(in) :: tproj
    real(kind=real64), intent(out) :: alpha, dlalpha
    real(kind=real64) :: df, dx, Lm, Lp, L0, alphamin, xmin
    integer(kind=int64) :: i

    call tstart("Lalpha_min")
    !backtracking
    dx = 1e-5
    alpha = 0.0D0
    call pgaxpy2(m, n, 1.0D0, x, alpha, s, xtmp)
    if (tproj) call doproj2(Xtmp, xmin)
    call L(Xtmp, L0)

    call pgaxpy2(m, n, 1.0D0, x, alpha + dx, s, xtmp)
    if (tproj) call doproj2(Xtmp, xmin)
    call L(Xtmp, Lp)
    df = (Lp - L0)/(dx)
    dlalpha = df
    alpha = 1.0D0

    lm = huge(lm)
    alphamin = 0

    do i = 1, 100
      call pgaxpy2(m, n, 1.0D0, x, alpha, s, xtmp)
      call doproj2(Xtmp, xmin)
      call L(Xtmp, Lp)
      if (lp .lt. lm .and. xmin .ge. -1e-12) then
        lm = lp
        alphamin = alpha
      end if

      if (rank .eq. 0) print *, "line search", i, "L=", Lp, "alpha=", alpha, "df(0)=", df, "min(X)=", xmin
      if (Lp .gt. L0 + 0.1D0*alpha*df .or. xmin .lt. -1e-12) then
        alpha = alpha*0.75D0
      else
        exit
      end if
    end do
    alpha = alphamin
    call tstop("Lalpha_min")
  end subroutine

  subroutine Lalpha_min_auglag(X, S, tproj, alpha, dlalpha)
    implicit none
    real(kind=real64), intent(inout), dimension(:, :) :: X, S
    logical, intent(in) :: tproj
    real(kind=real64), intent(out) :: alpha, dlalpha
    real(kind=real64) :: df, dx, Lm, Lp, L0, alphamin
    integer(kind=int64) :: i

    call tstart("Lalpha_min")
    !backtracking
    dx = 1e-5
    alpha = 0.0D0
    call pgaxpy2(m, n, 1.0D0, x, alpha, s, xtmp)
    if (tproj) call doproj(Xtmp)
    call L(Xtmp, L0)

    call pgaxpy2(m, n, 1.0D0, x, alpha + dx, s, xtmp)
    if (tproj) call doproj(Xtmp)
    call L(Xtmp, Lp)
    df = (Lp - L0)/(dx)
    dlalpha = df
    alpha = 1.0D0

    lm = huge(lm)
    alphamin = 0

    do i = 1, 100
      call pgaxpy2(m, n, 1.0D0, x, alpha, s, xtmp)
      if (tproj) call doproj(Xtmp)
      call L(Xtmp, Lp)
      if (lp .lt. lm) then
        lm = lp
        alphamin = alpha
      end if
      if (rank .eq. 0) print *, "line search", i, "L=", Lp, "alpha=", alpha, "df(0)=", df
      if (Lp .gt. L0 + 0.1D0*alpha*df) then
        alpha = alpha*0.75D0
      else
        exit
      end if
    end do
    alpha = alphamin
    call tstop("Lalpha_min")
  end subroutine

  subroutine paxpy(m, n, alpha, x, y)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :) :: x, y
    real(kind=real64), intent(in) :: alpha
    integer(kind=int64) :: i
    call tstart("paxpy")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      y(:, i) = y(:, i) + alpha*x(:, i)
    end do
    !$OMP end parallel DO
    call tstop("paxpy")
  end subroutine

  subroutine pgaxpy(m, n, alpha, x, beta, y)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :) :: x, y
    real(kind=real64), intent(in) :: alpha, beta
    integer(kind=int64) :: i
    call tstart("gpaxpy")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      y(:, i) = beta*y(:, i) + alpha*x(:, i)
    end do
    !$OMP end parallel DO
    call tstop("gpaxpy")
  end subroutine

  subroutine pgaxpy2(m, n, alpha, x, beta, y, z)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :) :: x, y, z
    real(kind=real64), intent(in) :: alpha, beta
    integer(kind=int64) :: i
    call tstart("pgpaxpy2")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      z(:, i) = beta*y(:, i) + alpha*x(:, i)
    end do
    !$OMP end parallel DO
    call tstop("pgpaxpy2")
  end subroutine

  subroutine pset(m, n, alpha, x)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(in) :: alpha
    real(kind=real64), intent(inout), dimension(:, :) :: x
    integer(kind=int64) :: i
    call tstart("pset")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      x(:, i) = alpha
    end do
    !$OMP end parallel DO
    call tstop("pset")
  end subroutine

  subroutine pcopy(m, n, x, y)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :) :: x, y
    integer(kind=int64) :: i
    call tstart("pcopy")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      y(:, i) = x(:, i)
    end do
    !$OMP end parallel DO
    call tstop("pcopy")
  end subroutine

  subroutine pscale(m, n, alpha, x)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :) :: x
    real(kind=real64), intent(in) :: alpha
    integer(kind=int64) :: i
    call tstart("pscale")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      x(:, i) = x(:, i)*alpha
    end do
    !$OMP end parallel DO
    call tstop("pscale")
  end subroutine

  subroutine pcopy_scale(m, n, alpha, x, y)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :) :: x, y
    real(kind=real64), intent(in) :: alpha
    integer(kind=int64) :: i
    call tstart("pcopy_scale")
    !$OMP parallel DO schedule(static,1)
    do i = 1, N
      y(:, i) = x(:, i)*alpha
    end do
    !$OMP end parallel DO
    call tstop("pcopy_scale")
  end subroutine

  subroutine pdot(m, n, x, y, val)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :), pointer :: x, y
    real(kind=real64), intent(out) :: val
    real(kind=real64), external :: ddot
    integer(kind=int64) :: i
#ifdef MPI
    integer(kind=int32) :: ierr
#endif
    call tstart("pdot")
    call tstart("pdot_red")
    val = 0
    !$OMP parallel DO schedule(static,1) reduction(+:val)
    do i = 1, N
      val = val + sum(x(:, i)*y(:, i))
    end do
    !$OMP end parallel DO
    call tstop("pdot_red")
#ifdef MPI
    call tstart("pdot_comm")
    call MPI_ALLREDUCE(MPI_IN_PLACE, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call tstop("pdot_comm")
#endif
    call tstop("pdot")
  end subroutine

  subroutine pmin(m, n, x, val)
    implicit none
    integer(kind=int64), intent(in) :: n, m
    real(kind=real64), intent(inout), dimension(:, :), pointer :: x
    real(kind=real64), intent(out) :: val
    integer(kind=int64) :: i
#ifdef MPI
    integer(kind=int32) :: ierr
#endif
    call tstart("pmin")
    call tstart("pmin_red")
    val = 0
    !$OMP parallel DO schedule(static,1) reduction(min:val)
    do i = 1, N
      val = minval(X(:, i))
    end do
    !$OMP end parallel DO
    call tstop("pmin_red")
#ifdef MPI
    call tstart("pmin_comm")
    call MPI_ALLREDUCE(MPI_IN_PLACE, val, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
    call tstop("pmin_comm")
#endif
    call tstop("pmin")
  end subroutine

  recursive subroutine butterfly(N, rank, offset, level, send, recv)
    implicit none
    integer(kind=int32), intent(in) :: N
    integer(kind=int32), intent(in) :: offset, level, rank
    integer(kind=int64), intent(inout), pointer, dimension(:, :) :: send, recv
    integer(kind=int32) :: N1, N2, a1, a2
    integer(kind=int32) :: i

    if (N .eq. 1) return

    !halves, N1<=N2
    N1 = N/2
    N2 = N - N1
    !print*,N,offset,level,N1,N2

    !everyone from half 1 needs to send to half 2
    do i = 1, N1
      a1 = i + offset
      a2 = N1 + i + offset
      if (a1 .ne. a2) then
        !print*,"send1",level,a1,a2
        if (a1 .eq. rank + 1) then
          if (send(level, 1) .eq. -1) then
            send(level, 1) = a2 - 1
          else
            send(level, 2) = a2 - 1
          end if
        end if
        if (a2 .eq. rank + 1) then
          if (recv(level, 1) .eq. -1) then
            recv(level, 1) = a1 - 1
          else
            recv(level, 2) = a1 - 1
          end if
        end if
      end if
    end do
    call butterfly(N1, rank, offset, level + 1, send, recv)

    !everyone from half 2 needs to send to half 1
    do i = N1 + 1, N
      a1 = i + offset
      a2 = mod(i + N1 + 1, N1) + 1 + offset
      if (a1 .ne. a2) then
        !print*,"send2",level,a1,a2
        if (a1 .eq. rank + 1) then
          if (send(level, 1) .eq. -1) then
            send(level, 1) = a2 - 1
          else
            send(level, 2) = a2 - 1
          end if
        end if
        if (a2 .eq. rank + 1) then
          if (recv(level, 1) .eq. -1) then
            recv(level, 1) = a1 - 1
          else
            recv(level, 2) = a1 - 1
          end if
        end if
      end if
    end do
    call butterfly(N2, rank, offset + N1, level + 1, send, recv)

  end subroutine

  subroutine redO(X, allreduce, O)
    implicit none
    real(kind=real64), pointer, intent(in), dimension(:, :) :: X
    logical, intent(in)  :: allreduce
    real(kind=real64), pointer, intent(in), dimension(:, :) :: O
    integer(kind=int64) :: i, j
#ifdef MPI
    TYPE(MPI_Status) :: stat
    logical :: found
    integer(kind=int32) :: ierr, rc
    real(kind=real64), allocatable, dimension(:) :: tmpf
    integer(kind=int64), pointer, dimension(:, :) :: send, recv
    integer(kind=int64) :: iop
    integer(kind=int32) :: k
    integer(kind=int32) :: lmax, level, lm, rcm, rankl
#endif

    call tstart("redO_O")
    !$OMP parallel DO private(j) schedule(static,1)
    do i = 1, D
      if (Flocalrow0(i) .ne. -1) then
        do J = 1, N
          O(j, i) = sum(Frows(i)%f*X(Flocalrow0(i):Flocalrow1(i), j))
        end do
      else
        O(:, i) = 0.0D0
      end if
    end do
    !$OMP end parallel DO
    call tstop("redO_O")
#ifdef MPI
    call tstart("redO_comm")
    !reduce O where necessary
    rc = 0
    if (.true.) then
      !send
      do i = 1, D
        if (block_red(i)%c .gt. 1) then
          found = .false.
          do k = 1, block_red(i)%c
            if (block_red(i)%f(k) .eq. rank) then
              found = .true.
              exit
            end if
          end do
          if (found) then
            if (allreduce) then
              do k = 1, block_red(i)%c
                if (block_red(i)%f(k) .ne. rank) then
                  rc = rc + 1
                  if (.not. requestsbuf(rc)%isalloc) then
                    allocate (requestsbuf(rc)%f(N))
                    requestsbuf(rc)%isalloc = .true.
                  end if
                  requestsbuf(rc)%f(:) = O(:, i)
                  call MPI_ISEND(requestsbuf(rc)%f, int(N, kind=int32), MPI_DOUBLE, &
                                 block_red(i)%f(k), int(i, kind=int32), MPI_COMM_WORLD, requests(rc), ierr)
!                  dsend(block_red(i)%f(k))=dsend(block_red(i)%f(k))+N*8
                end if
              end do
            else
              if (block_red(i)%f(1) .ne. rank) then
                rc = rc + 1
                if (.not. requestsbuf(rc)%isalloc) then
                  allocate (requestsbuf(rc)%f(N))
                  requestsbuf(rc)%isalloc = .true.
                end if
                requestsbuf(rc)%f = O(:, i)
                !              print*,rank,"MPI_ISEND",i,rank,block_red(i)%f(1)
                call MPI_ISEND(requestsbuf(rc)%f, int(N, kind=int32), MPI_DOUBLE, &
                               block_red(i)%f(1), int(i, kind=int32), MPI_COMM_WORLD, requests(rc), ierr)
!                dsend(block_red(i)%f(1))=dsend(block_red(i)%f(1))+N*8
              end if
            end if
          end if
        end if
      end do
      !receive
      allocate (tmpf(N))
      do i = 1, D
        if (block_red(i)%c .gt. 1) then
          found = .false.
          do k = 1, block_red(i)%c
            if (block_red(i)%f(k) .eq. rank) then
              found = .true.
              exit
            end if
          end do
          if (found) then
            if (allreduce) then
              do k = 1, block_red(i)%c
                if (block_red(i)%f(k) .ne. rank) then
!                  print*,rank,"MPI_RECV",i,block_red(i)%f(k),rank
                  call MPI_RECV(tmpf, int(N, kind=int32), MPI_DOUBLE, block_red(i)%f(k), &
                                int(i, kind=int32), MPI_COMM_WORLD, stat, ierr)
                  O(:, i) = O(:, i) + tmpf
!                  drecv(block_red(i)%f(k))=drecv(block_red(i)%f(k))+N*8
                end if
              end do
            else
              if (block_red(i)%f(1) .eq. rank) then
                do k = 1, block_red(i)%c
                  if (block_red(i)%f(k) .ne. rank) then
!                    print*,rank,"MPI_RECV",i,block_red(i)%f(k),rank
                    call MPI_RECV(tmpf, int(N, kind=int32), MPI_DOUBLE, block_red(i)%f(k), &
                                  int(i, kind=int32), MPI_COMM_WORLD, stat, ierr)
                    O(:, i) = O(:, i) + tmpf
!                    drecv(block_red(i)%f(1))=drecv(block_red(i)%f(1))+N*8
                  end if
                end do
              end if
            end if
          end if
        end if
      end do
      call mpi_waitall(rc, requests, MPI_STATUSES_IGNORE, ierr)
      deallocate (tmpf)
    else
      !butterfly communication for reduction
      allocate (tmpf(N))
      lm = 0
      do i = 1, D
        if (block_red(i)%c .le. 1) cycle
        lmax = bit_size(block_red(i)%c - 1) - leadz(block_red(i)%c - 1)
        lm = max(lm, lmax)
      end do

      allocate (send(lm, 2))
      allocate (recv(lm, 2))
      rcm = 0
      do level = 1, lm
        rc = 0
        do i = 1, D
          !check if this line needs a reduction
          if (block_red(i)%c .eq. 1) cycle

          !perform butterfly reduction
          lmax = bit_size(block_red(i)%c - 1) - leadz(block_red(i)%c - 1)
          if (level .gt. lmax) cycle

          !check if this rank is involved in the reduction
          found = .false.
          do k = 1, block_red(i)%c
            if (block_red(i)%f(k) .eq. rank) then
              found = .true.
              rankl = k - 1
              exit
            end if
          end do
          if (.not. found) cycle

          send(:, :) = -1
          recv(:, :) = -1
          call butterfly(block_red(i)%c, rankl, 0, 1, send, recv)

          do iop = 1, 2
            if (send(level, iop) .eq. -1) cycle
            rc = rc + 1
            if (.not. requestsbuf(rc)%isalloc) then
              allocate (requestsbuf(rc)%f(N))
              requestsbuf(rc)%isalloc = .true.
            end if
            requestsbuf(rc)%f(:) = O(:, i)
            call MPI_ISEND(requestsbuf(rc)%f, int(N, kind=int32), MPI_DOUBLE, &
                           block_red(i)%f(send(level, iop) + 1), int(i, kind=int32), MPI_COMM_WORLD, requests(rc), ierr)
          end do
        end do

        do i = 1, D
          !check if this line needs a reduction
          if (block_red(i)%c .eq. 1) cycle

          !perform butterfly reduction
          lmax = bit_size(block_red(i)%c - 1) - leadz(block_red(i)%c - 1)
          if (level .gt. lmax) cycle

          !check if this rank is involved in the reduction
          found = .false.
          do k = 1, block_red(i)%c
            if (block_red(i)%f(k) .eq. rank) then
              found = .true.
              rankl = k - 1
              exit
            end if
          end do
          if (.not. found) cycle

          send(:, :) = -1
          recv(:, :) = -1
          call butterfly(block_red(i)%c, rankl, 0, 1, send, recv)
          !recv
          do iop = 1, 2
            if (recv(level, iop) .eq. -1) cycle
            call MPI_RECV(tmpf, int(N, kind=int32), MPI_DOUBLE, &
                          block_red(i)%f(recv(level, iop) + 1), int(i, kind=int32), MPI_COMM_WORLD, stat, ierr)
            O(:, i) = O(:, i) + tmpf(:)
          end do
        end do
        rcm = max(rcm, rc)
      end do
      call mpi_waitall(rc, requests, MPI_STATUSES_IGNORE, ierr)
      deallocate (send)
      deallocate (recv)
      deallocate (tmpf)
    end if
    if (tbarrier) call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call tstop("redO_comm")
#endif
  end subroutine

  subroutine L(X, val)
    implicit none
    real(kind=real64), pointer, intent(in), dimension(:, :) :: X
    real(kind=real64), intent(out) :: val
    integer(kind=int64) :: i, j
#ifdef MPI
    integer(kind=int32) :: ierr
#endif

    call tstart("L")
    call tstart("L_prep")
    !reconstruct max column
    if (.not. tfullspace) then
      !$OMP parallel DO schedule(static,1)
      do j = 1, localrows
        X(j, maxr(j)) = 1.0D0 - (sum(X(j, :)) - X(j, maxr(j)))
      end do
      !$OMP end parallel DO
    end if
    call tstop("L_prep")
    call tstart("L_O")
    call redO(X, .true., O)
    call tstop("L_O")

    call tstart("L_red")
    call tstart("L_red_calc")
    val = 0.0D0
    !$OMP parallel DO reduction(+:val) schedule(static,1)
    do i = 1, D
      if (block_red(i)%c .gt. 0) then
        if (block_red(i)%f(1) .eq. rank) then
          val = val + sum((O(:, i) - P(:, i))**2)
        end if
      end if
    end do
    !$OMP end parallel DO
    call tstop("L_red_calc")
    if (gam .ne. 0) then
      call tstart("L_gam")
      !$OMP parallel DO reduction(+:val) schedule(static,1) private(j)
      do i = 1, N
        do j = 1, localrows - 1
          val = val + gam*(X(j, i) - X(j + 1, i))**2
        end do
      end do
      !$OMP end parallel DO
      call tstop("L_gam")
    end if

#ifdef MPI
    call tstart("L_red_comm")
    call MPI_ALLREDUCE(MPI_IN_PLACE, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call tstop("L_red_comm")
#endif
    call tstop("L_red")
    call tstop("L")
  end subroutine

  subroutine dL(X, dLdX)
    implicit none
    real(kind=real64), pointer, intent(in), dimension(:, :) :: X
    real(kind=real64), pointer, intent(inout), dimension(:, :) :: dLdX
    integer(kind=int64) :: i, j

    call tstart("dL")
    call tstart("dL_prep")
    !reconstruct max column
    if (.not. tfullspace) then
      !$OMP parallel DO schedule(static,1)
      do j = 1, localrows
        X(j, maxr(j)) = 1.0D0 - (sum(X(j, :)) - X(j, maxr(j)))
      end do
      !$OMP end parallel DO
    end if
    call pset(M, N, 0.0D0, dldx)
    call tstop("dL_prep")

    call tstart("dL_O")
    call redO(X, .true., O)
    call tstop("dL_O")
    call tstart("dL_dLdX")
    !$OMP parallel DO private(i) schedule(static,1)
    do j = 1, N
      do i = 1, D
        if (Flocalrow0(i) .ne. -1) then
          dLdX(Flocalrow0(i):Flocalrow1(i), j) = dLdX(Flocalrow0(i):Flocalrow1(i), j) + 2.0D0*(O(j, i) - P(j, i))*Frows(i)%f
        end if
      end do
    end do
    !$OMP end parallel DO
    call tstop("dL_dLdX")
    if (gam .ne. 0) then
      call tstart("dL_gam")
      !$OMP parallel DO schedule(static,1) private(j)
      do i = 1, N
        do j = 1, localrows
          if (j .lt. localrows) then
            dLdX(j, i) = dLdx(j, i) + 2.0D0*gam*(X(j, i) - X(j + 1, i))
          end if
          if (j .gt. 1) then
            dLdX(j, i) = dLdx(j, i) - 2.0D0*gam*(X(j - 1, i) - X(j, i))
          end if
        end do
      end do
      !$OMP end parallel DO
      call tstop("dL_gam")
    end if
    if (.not. tfullspace) then
      !$OMP parallel DO schedule(static,1)
      do j = 1, localrows
        dLdX(j, :) = dldx(j, :) - dLdX(j, maxr(j))
      end do
      !$OMP end parallel DO
    end if
    call tstop("dL")
  end subroutine

  subroutine conv_test(X, DLDX, tprint, conv)
    implicit none
    real(kind=real64), pointer, intent(in), dimension(:, :) :: X
    real(kind=real64), pointer, intent(inout), dimension(:, :) :: dLdX
    integer(kind=int64) :: i, j
    logical, intent(in) :: tprint
    real(kind=real64), intent(out) :: conv
    real(kind=real64) :: lambda

    !kkt first order conditions for X
    !reconstruct max column
    if (.not. tfullspace) then
      !$OMP parallel DO schedule(static,1)
      do j = 1, localrows
        X(j, maxr(j)) = 1.0D0 - (sum(X(j, :)) - X(j, maxr(j)))
      end do
      !$OMP end parallel DO
    end if
    tfullspace = .true.
    call dl(X, dldX)
    conv = 0
    !$OMP parallel DO private(lambda,j) reduction(+:conv)
    do i = 1, localrows
      lambda = max(0.0D0, -minval(dldx(i, :)))
      conv = conv + sum(((dldX(i, :) + lambda)*X(i, :))**2)
      if (tprint) print *, "lambda", i, lambda, minval(dldX(i, :) + lambda)
    end do
    !$OMP end parallel DO
#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, conv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    conv = sqrt(conv/(real(N, kind=real64)*real(M, kind=real64)))
  end subroutine

  subroutine read_unformatted_binary_from_python(M, N, filename, A)
    implicit none
    integer(kind=int64), intent(in)  :: M, N
    real(kind=real64), pointer, intent(in), dimension(:, :) :: A
    character(len=*), intent(in) :: filename

    open (40, file=filename, status='old', access='stream', form='unformatted')
    read (40) A
    close (40)
  end subroutine

  subroutine read_sparse_from_python(L, f1, f2, f3, Ai, Aj, Ad)
    implicit none
    integer(kind=int64), intent(in)  :: L
    integer(kind=int64), pointer, dimension(:) :: Ai, Aj
    real(kind=real64), pointer, dimension(:) :: Ad
    character(len=*), intent(in) :: f1, f2, f3
    integer(kind=int32), pointer, dimension(:) :: Ai2, Aj2

    allocate (Ai(L))
    allocate (Aj(L))
    allocate (Ai2(L))
    allocate (Aj2(L))
    allocate (Ad(L))

    open (40, file=f1, status='old', access='stream', form='unformatted')
    read (40) Ai2
    close (40)
    Ai2 = Ai2 + 1
    ai(:) = ai2(:)
    open (40, file=f2, status='old', access='stream', form='unformatted')
    read (40) Aj2
    close (40)
    Aj2 = Aj2 + 1
    aj(:) = aj2(:)
    open (40, file=f3, status='old', access='stream', form='unformatted')
    read (40) Ad
    close (40)
  end subroutine

  subroutine set_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer(kind=int32) :: n

    call random_seed(size=n)
    allocate (seed(n))
    seed(:) = 23764
    call random_seed(put=seed)
  end subroutine

  subroutine doproj(X)
    implicit none
    real(kind=real64), pointer, intent(inout), dimension(:, :) :: X
    integer(kind=int64) :: i
    real(kind=real64), dimension(:), allocatable :: x2, x3
    integer(kind=int32) :: n2
    call tstart("doproj")
    allocate (x2(N))
    allocate (x3(N))
    n2=int(N,kind=int32)
    !$OMP parallel DO private(x2,x3) schedule(static,1)
    do i = 1, localrows
      x3 = X(i, :)
      call simplexproj_Condat(x3, x2, n2)
      X(i, :) = x2
    end do
    !$OMP end parallel DO
    call tstop("doproj")
    deallocate (x2)
    deallocate (x3)
  end subroutine

  subroutine doproj2(X, xmin)
    implicit none
    real(kind=real64), pointer, intent(inout), dimension(:, :) :: X
    real(kind=real64), intent(out) :: xmin
    integer(kind=int64) :: i, j
    call tstart("doproj2")
    xmin = 0
    !$OMP parallel DO private(j) schedule(static,1) reduction(min:xmin)
    do i = 1, localrows
      do j = 1, N
        if (X(i, j) .lt. 0) then
          X(i, j) = 0.0D0
        end if
      end do
      X(i, maxr(i)) = 1.0D0 - (sum(x(i, :)) - X(i, maxr(i)))
      xmin = min(xmin, X(i, maxr(i)))
    end do
    !$OMP end parallel DO
#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
    call tstop("doproj2")
  end subroutine

  function rowdist(i) result(j)
    integer(kind=int64), intent(in) :: i
    integer(kind=int64) :: j
    if (i .le. MB1*NP1) then
      j = min((i - 1)/MB1, num_proc - 1)
    else
      j = NP1 + ((i - MB1*NP1) - 1)/MB2
    end if
  end function

  function map_global_to_local(i) result(j)
    integer(kind=int64), intent(in) :: i
    integer(kind=int64) :: j
    if (rowdist(i) .eq. rank) then
      if (i .le. MB1*NP1) then
        j = mod(i - 1, MB1) + 1
      else
        j = mod(i - MB1*NP1 - 1, MB2) + 1
      end if
    else
      j = 0
    end if
  end function
end module

program pqdts
  use basics
  implicit none
  integer(kind=int64) :: i, j, k, i0, j0, i1, il0, il1, ic, dend, dstart
  CHARACTER(len=32) :: arg
  real(kind=real64), pointer, dimension(:, :) :: X, dLdX, z, bd
  integer(kind=int64), pointer, dimension(:) :: tmp, tmp2
  integer(kind=int32) :: ntmp
  real(kind=real64) :: val, t0, t1, t2, l0, val2, ts, svar
  real(kind=real64) :: alpha, beta, b1, b2, norm_der, v2, eps
  real(kind=real64) :: pos, fmean
  real(kind=real64) :: conv, xmin
  logical :: tproj = .true.
  logical, parameter :: testder = .false.
  logical, parameter :: tbench = .false.
  real(kind=real64), parameter :: dxp = 1e-7
  real(kind=real64) :: tactivetol = 1e-8
  real(kind=real64) :: conv_thres = 1e-6
  logical :: found
  integer(kind=int32) :: iiter, maxiiter, oiter, iter2, smo_dist, smo_c
  integer(kind=int64) :: mean
  Character(len=255) :: outname
  logical :: file_exists
  integer(kind=int32) :: tactive = 0
  logical :: precond = .false.
  real(kind=real64) :: Fsize = 0, dlalpha, smo_fact
  integer(kind=int32) :: output, state, start_stage, rank2, read_input
  integer(kind=int64) :: nact
  integer(kind=int64) :: M2, N2, D2, localrows2
  Character(len=255) :: str

#ifdef MPI
  tbarrier = 1 .eq. 1
#endif
  tfullspace = .true.

  dlalpha = huge(dlalpha)

  allocate (timingt(50))
  timingt(:) = 0.0D0
  allocate (timingsum(50))
  timingsum(:) = 0.0D0
  allocate (timingc(50))
  timingc(:) = 0
  allocate (timingn(50))
  timingid = 0

  call tstart("prep")

  t0 = omp_get_wtime()
  call set_seed()

  CALL getarg(1, arg)
  read (arg, *) M
  CALL getarg(2, arg)
  read (arg, *) N
  CALL getarg(3, arg)
  read (arg, *) D
  CALL getarg(4, arg)
  read (arg, *) Fl
  CALL getarg(5, arg)
  read (arg, *) Pl
  CALL getarg(6, arg)
  read (arg, *) mode
  CALL getarg(7, arg)
  read (arg, *) maxiiter
  CALL getarg(8, arg)
  read (arg, *) output
  CALL getarg(9, arg)
  read (arg, *) gam
  CALL getarg(10, arg)
  read (arg, *) conv_thres
  CALL getarg(11, arg)
  read (arg, *) start_stage
  CALL getarg(12, arg)
  read (arg, *) read_input
  CALL getarg(13, arg)
  read (arg, *) smo_fact

#ifdef MPI
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
#else
  rank = 0
  num_proc = 1
#endif
  if (abs(gam) .ne. 0 .and. num_proc .gt. 1) then
    print *, "WARNING: smoothing currently only correct for single-rank because boundaries are not accounted for."
  end if

  allocate (dsend(num_proc))
  dsend(:) = 0
  allocate (drecv(num_proc))
  drecv(:) = 1

  if (rank .eq. 0) print *, "num_proc", num_proc

  !determine distribution of rows

  MB = M/num_proc
  MB1 = MB
  MB2 = MB1 + 1
  NP2 = M - num_proc*MB
  NP1 = num_proc - NP2
  if (rank .eq. 0) print *, "M=", M, "N=", N, "D=", D, "MB=", MB, "MB1=", MB1, "MB2=", MB2, "NP1=", NP1, "NP2=", NP2
  if (rank .eq. 0) print *, "mode=", mode, "maxiiter=", maxiiter, "output=", output
  if (rank .eq. 0) print *, "gam=", gam, "smo_fact=", smo_fact
  if (rank .eq. 0) print *, "conv_thres=", conv_thres

  if (rank + 1 .le. NP1) then
    localrows0 = MB1*(rank) + 1
    localrows1 = localrows0 + MB1 - 1
    localrows = MB1
  else
    localrows0 = MB1*NP1 + (rank - NP1)*MB2 + 1
    localrows1 = localrows0 + MB2 - 1
    localrows = MB2
  end if

  print *, rank, "localrows", localrows0, localrows1, localrows, rowdist(localrows0), rowdist(localrows1)

  if (Fl .ge. 0) then
    !FIXME only read F partially
    call read_sparse_from_python(Fl, "data/F_row.bin", "data/F_col.bin", "data/F_data.bin", Fi, Fj, Fd)
    print *, rank, "sum(F)=", sum(Fd), minval(Fi), maxval(Fi), minval(Fj), maxval(Fj)

    allocate (Frow0(D))
    allocate (Frow1(D))
    allocate (Flocalrow0(D))
    allocate (Flocalrow1(D))
    allocate (Frows(D))
    !$OMP parallel DO private(i0,i1,il0,il1,j) schedule(static,1)
    do i = 1, d
      i0 = huge(i0)
      i1 = 0
      il0 = huge(i0)
      il1 = 0
      do j = 1, Fl
        if (Fi(j) .eq. i) then
          i0 = min(Fj(j), i0)
          i1 = max(Fj(j), i1)
          if (map_global_to_local(Fj(j)) .ne. 0) then
            il0 = min(Fj(j), il0)
            il1 = max(Fj(j), il1)
          end if
        end if
      end do
      !if(rank.eq.0)print*,"Rows in F",rank,i,i0,i1,i1-i0+1,il0,il1,il1-il0+1
      Frow0(i) = i0
      Frow1(i) = i1
      if (il0 .ne. huge(i0)) then
        Flocalrow0(i) = map_global_to_local(il0)
      else
        Flocalrow0(i) = -1
      end if
      if (il1 .ne. 0) then
        Flocalrow1(i) = map_global_to_local(il1)
      else
        Flocalrow1(i) = -1
      end if
      allocate (Frows(i)%f(il1 - il0 + 1))
      Frows(i)%f(:) = 0
      do j = 1, Fl
        if (Fi(j) .eq. i .and. map_global_to_local(Fj(j)) .ne. 0) then
          Frows(i)%f(Fj(j) - il0 + 1) = Fd(j)
        end if
      end do
    end do
    deallocate (Fd, Fi, Fj)
  else
    !build rows of F internally
    !M = (D-1)**2
    !for i in range(D):
    !  mean = (D-i-1)**2
    !  if (mean < 50):
    !    row = sp.stats.poisson.pmf(range(M), mean)
    !  else:
    !    row = sp.stats.norm.pdf(range(M), mean, np.sqrt(mean))
    if (M .ne. (D - 1)**2) then
      print *, "M is not (D-1)**2"
      stop
    end if
    allocate (Frow0(D))
    allocate (Frow1(D))
    Frow0(:) = 0
    Frow1(:) = 0

    dstart = 1 + rank*(D/num_proc + 1)
    dend = min(D, dstart + (D/num_proc + 1) - 1)
    !$OMP parallel DO private(mean,i0,i1,il0,il1,fmean,j,pos,svar) schedule(static,1)
    do i = dstart, dend
      mean = (D - i)**2
      i0 = huge(i0)
      i1 = 0
      il0 = huge(i0)
      il1 = 0
      do j = max(1, mean), M
        pos = real(j - 1, kind=real64)
        fmean = real(mean, kind=real64)
        call fprob(pos, fmean, svar)
        if (svar .lt. 1d-10) exit
        i0 = min(j, i0)
        i1 = max(j, i1)
      end do
      do j = max(1, mean), 1, -1
        pos = real(j - 1, kind=real64)
        fmean = real(mean, kind=real64)
        call fprob(pos, fmean, svar)
        if (svar .lt. 1d-10) exit
        i0 = min(j, i0)
        i1 = max(j, i1)
      end do
      Frow0(i) = i0
      Frow1(i) = i1
    end do
    !$OMP end parallel DO
#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, Frow0, int(D, kind=int32), MPI_LONG, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, Frow1, int(D, kind=int32), MPI_LONG, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

    allocate (Flocalrow0(D))
    allocate (Flocalrow1(D))
    allocate (Frows(D))
    fsize = 0

    !$OMP parallel DO private(mean,i0,i1,il0,il1,fmean,j,pos,svar) reduction(+:fsize) schedule(static,1)
    do i = 1, D
      mean = (D - i)**2
      il0 = huge(i0)
      il1 = 0

      if (localrows0 .le. Frow0(i) .and. localrows1 .ge. Frow0(i)) then
        il0 = Frow0(i)
        il1 = min(localrows1, Frow1(i))
      end if
      if (localrows0 .le. Frow1(i) .and. localrows1 .ge. Frow1(i)) then
        il0 = max(localrows0, Frow0(i))
        il1 = Frow1(i)
      end if
      if (localrows0 .ge. Frow0(i) .and. localrows1 .le. Frow1(i)) then
        il0 = localrows0
        il1 = localrows1
      end if

      if (il0 .ne. huge(i0)) then
        Flocalrow0(i) = map_global_to_local(il0)
      else
        Flocalrow0(i) = -1
      end if
      if (il1 .ne. 0) then
        Flocalrow1(i) = map_global_to_local(il1)
      else
        Flocalrow1(i) = -1
      end if
      if (il0 .ne. 0) then
        allocate (Frows(i)%f(il1 - il0 + 1))
        Fsize = Fsize + (il1 - il0 + 1)*8
        Frows(i)%f(:) = 0
        do j = il0, il1
          pos = real(j - 1, kind=real64)
          fmean = real(mean, kind=real64)
          call fprob(pos, fmean, svar)
          Frows(i)%f(j - il0 + 1) = svar
        end do
      end if
    end do
    !$OMP end parallel DO
    print *, "Fsize", rank, fsize
  end if
  !do i=1,D
  !  write(44,*)"F",Flocalrow0(i),Flocalrow1(i),sum(Frows(i)%f)
  !enddo
  !call flush(44)

  !determine which nodes have blocks that need to be reduced
  allocate (block_red(D))
  ntmp = 10
  allocate (tmp(ntmp))
  do j = 1, D
    if (Flocalrow0(j) .eq. -1) then
      cycle
    end if
    ic = 0
    do i = Frow0(j), Frow1(j)
      found = .false.
      do k = 1, ic
        if (tmp(k) .eq. rowdist(i)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        ic = ic + 1
        if (ic .gt. ntmp) then
          allocate (tmp2(ntmp))
          tmp2(:) = tmp(:)
          deallocate (tmp)
          allocate (tmp(ntmp + 10))
          tmp(1:ntmp) = tmp2(1:ntmp)
          ntmp = ntmp + 10
        end if
        tmp(ic) = rowdist(i)
      end if
    end do
    allocate (block_red(j)%f(ic))
    block_red(j)%f(:) = int(tmp(1:ic), kind=int32)

    found = .false.
    do k = 1, ic
      if (tmp(k) .eq. rank) then
        found = .True.
        exit
      end if
    end do
    do k = 1, ic
      if (ic .gt. 1 .and. tmp(k) .ne. rank .and. found) then
        blocks_red = blocks_red + 1
      end if
    end do
    block_red(j)%c = int(ic, kind=int32)
  end do
  deallocate (tmp)
  print *, rank, "blocks_red", blocks_red, ntmp

#ifdef MPI
  allocate (requests(blocks_red))
  allocate (requestsbuf(blocks_red))
#endif

  !build P
  allocate (P(N, D))
  if (Pl .ge. 0) then
    call read_sparse_from_python(Pl, "data/P_row.bin", "data/P_col.bin", "data/P_data.bin", Pi, Pj, Pd)
    print *, rank, "sum(P)=", sum(Pd), minval(Pi), maxval(Pi), minval(Pj), maxval(Pj)
    P(:, :) = 0
    do i = 1, Pl
      P(Pj(i), Pi(i)) = Pd(i)
    end do
  else
    call random_number(P)
  end if

  allocate (X(localrows, N))
  call pset(M, N, 0.0D0, x)
  allocate (dLdX(localrows, N))
  call pset(M, N, 0.0D0, dldx)
  allocate (bd(localrows, N))
  call pset(M, N, 0.0D0, bd)
  !temporary
  allocate (O(N, D))
  allocate (Xtmp(localrows, N))
  call pset(M, N, 0.0D0, xtmp)
  allocate (sn(localrows, N))
  call pset(M, N, 0.0D0, sn)
  allocate (z(localrows, N))
  call pset(M, N, 0.0D0, z)

!  allocate (btmp(localrows, N))
!  call pset(M, N, 0.0D0, btmp)
!  allocate (btmp2(localrows, N))
!  call pset(M, N, 0.0D0, btmp2)
!  allocate (btmp3(localrows, N))
!  call pset(M, N, 0.0D0, btmp3)

  allocate (ik(localrows, N))
  allocate (maxr(localrows))

  do i = 1, localrows
    call random_number(X(i, :))
!    X(i, :) = 1.0D0
    X(i, :) = X(i, :)/sum(X(i, :))
  end do

  !$OMP parallel DO schedule(static,1)
  do j = 1, localrows
    maxr(j) = maxloc(X(j, :), 1)
  end do
  !$OMP end parallel DO

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  if (tbench) then
    if (rank .eq. 0) print *, "quick benchmark:"
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t1 = omp_get_wtime()
      call L(X, val)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t2 = omp_get_wtime()
      if (rank .eq. 0) print *, "t L", (t2 - t1), val
    end do

    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t1 = omp_get_wtime()
      call dL(X, dLdX)
      t2 = omp_get_wtime()
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t dL", t2 - t1, dLdX(1, 1)
    end do

    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pset(M, N, 1.0D0, dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pset", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pcopy(M, N, X, dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pcopy", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pscale(M, N, -1.0D0, dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pscale", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pcopy_scale(M, N, -1.0D0, X, dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pcopy_scale", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call paxpy(M, N, 1.0D0, X, dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t paxpy", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pgaxpy(M, N, 2.0D0, X, 0.5D0, dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pgaxpy", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pgaxpy2(M, N, 2.0D0, X, 0.5D0, dLdX, bd)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pgaxpy2", omp_get_wtime() - ts
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call pdot(M, N, X, Dldx, svar)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t pdot", omp_get_wtime() - ts, svar
    end do
    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ts = omp_get_wtime()
      call doproj(dLdX)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      if (rank .eq. 0) print *, "t doproj", omp_get_wtime() - ts
    end do

    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t1 = omp_get_wtime()
      call Bdiag(bd)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t2 = omp_get_wtime()
      if (rank .eq. 0) print *, "t Bdiag", (t2 - t1), sum(bd(1:4, 1:N))
    end do

    do i = 1, 10
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t1 = omp_get_wtime()
      call BF(dLdX, bd)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t2 = omp_get_wtime()
      if (rank .eq. 0) print *, "t BF", (t2 - t1), sum(bd(1:2, 1:N))
    end do
  end if

  !test derivatives
  if (testder) then
    allocate (Xp(localrows, N))
    call pset(M, N, 0.0D0, xp)
    if (num_proc .gt. 1) then
      print *, "Derivative test only for a single mpi-rank"
      stop
    end if
    do ic = 1, 2
      if (ic .eq. 1) tfullspace = .true.
      if (ic .eq. 2) tfullspace = .false.
      if (.not. tfullspace) then
        !$OMP parallel DO schedule(static,1)
        do j = 1, localrows
          X(j, maxr(j)) = 1.0D0 - (sum(X(j, :)) - X(j, maxr(j)))
        end do
        !$OMP end parallel DO
      end if

      if (rank .eq. 0) print *, "testing derivatives: tfullspace=", tfullspace
      call pset(M, N, 0.0D0, xp)
      call dL(X, dLdX)
      call L(X, L0)
      do i0 = -2, 3 !M
        do j0 = -2, 3 !N
          if (i0 .le. 0) then
            i = M + i0
          else
            i = i0
          end if
          if (j0 .le. 0) then
            j = N + j0
          else
            j = j0
          end if
          xp(:, :) = X
          if (map_global_to_local(i) .ne. 0) Xp(map_global_to_local(i), j) = Xp(map_global_to_local(i), j) + dxp
          call L(Xp, val)
          xp(:, :) = X
          if (map_global_to_local(i) .ne. 0) Xp(map_global_to_local(i), j) = Xp(map_global_to_local(i), j) - dxp
          call L(Xp, val2)
          if (rank .eq. 0) print *, "derivative test", i, j, (val - val2)/(2.0D0*dxp), dLdX(map_global_to_local(i), j), &
            (val - val2)/(2.0D0*dxp) - dLdX(map_global_to_local(i), j)
        end do
      end do

      if (rank .eq. 0) print *, "testing second derivatives times vector: tfullspace=", tfullspace
      call random_number(sn)

      call BF(sn, bd)
      call pset(M, N, 0.0D0, xp)
      call dl(X, dldx)
      xp = X + dxp*sn
      call dl(Xp, sn)

      do i0 = -2, 3 !M
        do j0 = -2, 3 !N
          if (i0 .le. 0) then
            i = M + i0
          else
            i = i0
          end if
          if (j0 .le. 0) then
            j = N + j0
          else
            j = j0
          end if
          if (rank .eq. 0) print *, "second derivative test", i, j, (sn(i, j) - dldx(i, j))/dxp, bd(i, j), &
            (sn(i, j) - dldx(i, j))/dxp - bd(i, j)
        end do
      end do

      if (rank .eq. 0) print *, "testing second derivatives diagonal: tfullspace=", tfullspace
      call pset(M, N, 0.0D0, xp)
      call bdiag(bd)
      do i0 = -2, 3 !M
        do j0 = -2, 3 !N
          if (i0 .le. 0) then
            i = M + i0
          else
            i = i0
          end if
          if (j0 .le. 0) then
            j = N + j0
          else
            j = j0
          end if
          call pset(M, N, 0.0D0, xp)
          xp(i, j) = 1.0D0
          call BF(xp, sn)
          if (rank .eq. 0) print *, "second derivative diagonal test", i, j, sn(i, j), bd(i, j), sn(i, j) - bd(i, j)
        end do
      end do
    end do
    deallocate (Xp)
  end if

  !initial
  do j = 1, localrows
    X(j, :) = 1.0/real(N, kind=real64)
    !call random_number(X(j, :))
    !X(j, :) = X(j, :)/sum(X(j, :))
  end do

  !read inital state if it exists
  oiter = start_stage
  INQUIRE (FILE=trim(str), EXIST=file_exists)
  write (outname, '(a,i6,a,i6,a)') "rank_", rank, "_oiter", oiter, ".dat"
  INQUIRE (FILE=outname, EXIST=file_exists)
  if (file_exists .and. read_input .eq. 1) then
    print *, "reading starting state from", trim(adjustl(outname))
    call tstart("input")
    open (42, file=adjustl(trim(outname)), status='old', FORM='unformatted')
    read (42) M2
    read (42) N2
    read (42) D2
    read (42) localrows2
    read (42) rank2
    if (M .ne. M2) then
      print *, "input file in this directory is for a different M", M, "!=", M2
      stop
    end if
    if (N .ne. N2) then
      print *, "input file in this directory is for a different N", N, "!=", N2
      stop
    end if
    if (D .ne. D2) then
      print *, "input file in this directory is for a different D", D, "!=", D2
      stop
    end if
    if (localrows .ne. localrows2) then
      print *, "input file in this directory is for a different localrows", localrows, "!=", localrows2
      stop
    end if
    if (rank .ne. rank2) then
      print *, "input file in this directory is for a different rank", rank, "!=", rank2
      stop
    end if
    read (42) X
    close (42)
    call tstop("input")
    !smoothen
    if (localrows .ne. M .and. smo_fact .gt. 0) then
      print *, "smoothing step is currently only available for single-rank calcuclations"
      stop
    end if
    if (smo_fact .gt. 0) then
      call pcopy(M, N, x, dldx)
      !$OMP parallel DO schedule(static,1) private(i,smo_c,smo_dist)
      do j = 1, M
        smo_dist = int(anint(j/smo_fact), kind=int32)
        if (j .le. 100) cycle
        do i = 1, N
          smo_c = 0
          dldx(j, i) = 0
          do k = max(1, j - smo_dist), min(j + smo_dist, M)
            dldx(j, i) = dldx(j, i) + x(k, i)
            smo_c = smo_c + 1
          end do
          dldx(j, i) = dldx(j, i)/(smo_c)
        end do
        dldx(j, :) = dldx(j, :)/sum(dldx(j, :))
      end do
      !$OMP end parallel DO

      call pcopy(M, N, dldx, x)
    end if
  end if

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  call L(X, L0)
  call dL(X, dLdX)
  call pdot(M, N, dLdX, dLdX, norm_der)

  iiter = 0
  oiter = 0
  if (rank .eq. 0) print *, "iter", oiter, iiter, "L=", L0, &
    "||der||=", sqrt(norm_der), "dt=", t2 - t1, "ttot=", omp_get_wtime() - t0

  call tstop("prep")
  call tprint(0)
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  call tstart("min")

  !minimization
  iter2 = 0

  do oiter = start_stage, 2
    dlalpha = -huge(dlalpha)
    do iiter = 1, maxiiter
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      t1 = omp_get_wtime()
      if (mode .eq. 1) then
        !projected gradient
        tactive = 0
        tfullspace = .false.
        if (.not. tfullspace) then
          !transform to representation with implicit sum consrtaint
          !$OMP parallel DO schedule(static,1)
          do j = 1, localrows
            maxr(j) = maxloc(X(j, :), 1)
            X(j, maxr(j)) = 1.0D0
          end do
          !$OMP end parallel DO
        else
          maxr(:) = -1
        end if
        call dL(X, dLdX)
        call pcopy_scale(M, N, -1.0D0, dLDx, sn)
        call Lalpha_min2(X, sn, tproj, alpha, dlalpha)

      else if (mode .eq. 2) then
        !truncated newton
        if ((abs(dlalpha) .lt. 1e-4 .or. dlalpha .gt. 0) .and. oiter .eq. 1) then
          if (rank .eq. 0) print *, "stopping inner iteration because dldlpha is small", dlalpha
          exit
        end if
        if (oiter .eq. 1) then
          tactive = 1
          precond = .true.
          tfullspace = .true.
        end if
        if (oiter .eq. 2) then
          tactive = 1
          precond = .true.
          tfullspace = .false.
        end if

        if (.not. tfullspace) then
          !transform to representation with implicit sum consrtaint
          !$OMP parallel DO schedule(static,1)
          do j = 1, localrows
            if (maxr(j) .ne. -1) then
              X(j, maxr(j)) = 1.0D0 - (sum(X(j, :)) - X(j, maxr(j)))
            end if
            maxr(j) = maxloc(X(j, :), 1)
            X(j, maxr(j)) = 1.0D0
          end do
          !$OMP end parallel DO
        else
          maxr(:) = -1
        end if

        call dL(X, dLdX)
        if (tactive .ne. 0) then
          call tstart("active_dLdX")
          nact = 0
          !$OMP parallel DO schedule(static,1) private(j) reduction(+:nact)
          do i = 1, N
            do j = 1, localrows
              if (X(j, i) .lt. tactivetol .and. dldx(j, i) .gt. 0) then
                dLdX(j, i) = 0
                ik(j, i) = 1
                nact = nact + 1
              else
                ik(j, i) = 0
              end if
            end do
          end do
          !$OMP end parallel DO
          if (rank .eq. 0) print *, "positivity constraints", nact, "active of", M*N, &
            real(nact, kind=real64)/real(M*N, kind=real64)
          call tstop("active_dLdX")
        end if

        call pdot(M, N, dLdX, dLdX, norm_der)
        beta = 0
        alpha = 0

        call pset(m, n, 0.0D0, sn)
        if (precond) then
          call pset(m, n, 0.0D0, z)
        end if

        call pscale(M, N, -1.0D0, dLdX)

        if (precond) then
          call Bdiag(bd)
          call tstart("precond1")
          !$OMP parallel DO schedule(static,1)
          do i = 1, N
            z(:, i) = dLdX(:, i)/sqrt(bd(:, i)) !dLdX=r
          end do
          !$OMP end parallel DO
          call tstop("precond1")
        end if
        eps = min(0.5, sqrt(sqrt(norm_der)))*sqrt(norm_der)

        if (precond) then
          call pcopy(M, N, z, xtmp)
        else
          call pcopy(M, N, dLdX, xtmp) !dLdX=r
        end if

        do j = 1, 50
          call BF(xtmp, bd)

          call tstart("CG")
          if (tactive .ne. 0) then
            call tstart("active")
            !$OMP parallel DO schedule(static,1) private(k)
            do i = 1, N
              do k = 1, localrows
                if (ik(k, i) .eq. 1) then
                  bd(k, i) = 0
                end if
              end do
            end do
            !$OMP end parallel DO
            call tstop("active")
          end if
          if (precond) then
            call pdot(m, n, dLdX, z, b1) !dLdX=r
          else
            call pdot(m, n, dLdX, dLdX, b1) !dLdX=r
          end if
          call pdot(m, n, xtmp, bd, v2)

          call paxpy(m, n, b1/v2, xtmp, sn)
          call paxpy(m, n, -b1/v2, bd, dLdX) !dLdX=r

          if (precond) then
            call tstart("precond2")
            call bdiag(bd)
            !$OMP parallel DO schedule(static,1)
            do i = 1, N
              z(:, i) = dLdX(:, i)/sqrt(bd(:, i)) !dLdX=r
            end do
            !$OMP end parallel DO
            call tstop("precond2")
          end if
          if (precond) then
            call pdot(m, n, dLdX, z, v2) !dLdX=r
            call pgaxpy(M, N, 1.0D0, z, v2/b1, xtmp)
            call pdot(m, n, dLdX, dLdX, b2) !dLdX=r
          else
            call pdot(m, n, dLdX, dLdX, v2) !dLdX=r
            call pgaxpy(M, N, 1.0D0, dLdX, v2/b1, xtmp) !dLdX=r
            !call pdot(m,n,dLdX,dLdX,b2) !dLdX=r
            b2 = v2
          end if

          if (rank .eq. 0) print *, "newton-pcg", j, "sqrt(b2)=", sqrt(b2), "eps=", eps, "b1=", b1, "v2=", v2
          call tstop("CG")

          if (sqrt(b2) .le. eps) then
            exit
          end if
        end do

        if (1 .eq. 1) then
          call dL(X, dLdX)
          call bdiag(bd)
          !$OMP parallel DO schedule(static,1) private(k)
          do i = 1, N
            do k = 1, localrows
              if (ik(k, i) .eq. 1 .and. maxr(k) .ne. i) then
                sn(k, i) = -dldx(k, i)/bd(k, i)
              end if
            end do
          end do
          !$OMP end parallel DO
        end if
        if (tfullspace) then
          call Lalpha_min_auglag(X, sn, tproj, alpha, dlalpha)
        else
          call Lalpha_min2(X, sn, tproj, alpha, dlalpha)
        end if
      end if
      if (alpha .ne. 0) then
        call paxpy(m, n, alpha, sn, X)
      end if

      xmin = 0
      if (tproj) then
        if (tfullspace) then
          call doproj(X)
        else
          call doproj2(X, xmin)
        end if
        call L(X, L0)
        call dL(X, dLdX)
        call pdot(M, N, dLdX, dLdX, norm_der)
      end if
      if (tproj) then
        call conv_test(X, DLDX, .false., conv)
      end if

      t2 = omp_get_wtime()
      if (rank .eq. 0) print *, "iter", oiter, iiter, "L=", L0, "max(abs(c))", xmin, "conv", conv, &
        "||der||=", sqrt(norm_der), "mode=", mode, "alpha=", alpha, "dlalpha=", dlalpha, &
        "dt=", t2 - t1, "ttot=", omp_get_wtime() - t0
      if (conv .lt. conv_thres) then
        exit
      end if

    end do
    call conv_test(X, DLDX, output.ne.0, conv)
    if (output .ge. 1) then
      call tstart("output")
      !export
      write (outname, '(a,i6,a,i6,a)') "rank_", rank, "_oiter", oiter, ".dat"
      print *, outname
      INQUIRE (FILE=outname, EXIST=file_exists)
      if (file_exists) then
        open (42, file=trim(outname), status='old', FORM='unformatted')
      else
        open (42, file=trim(outname), status='new', FORM='unformatted')
      end if
      write (42) M
      write (42) N
      write (42) D
      write (42) localrows
      write (42) rank
      write (42) X
      close (42)
      call tstop("output")
    end if
    call tprint(oiter)
  end do
  call tstop("min")
  call tprint(-1)

  !get memory usage
  open (43, file="/proc/self/status", status='old', iostat=state)
  write (str, *) "status_", rank, ".out"
  inquire (file=trim(adjustl(str)), exist=file_exists)
  if (file_exists) then
    open (44, file=trim(adjustl(str)), status="replace", action="write")
  else
    open (44, file=trim(adjustl(str)), status="new", action="write")
  end if
  DO
    read (43, *, iostat=state) outname, arg
    if (state /= 0) EXIT
    write (44, *) trim(adjustl(outname)), " ", trim(adjustl(arg))
  end do
  close (43)
  close (44)

#ifdef MPI
  call MPI_Finalize(ierror)
#endif
end program
