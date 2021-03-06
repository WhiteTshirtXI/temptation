subroutine poisson(p, ut, vt, dx, dt, beta, tol, ni, nj)
  implicit none
  integer :: ni, nj, i, j
  real(8), intent(in) :: beta, tol, dx, dt
  real(8), intent(in), dimension(ni+1, nj) :: ut
  real(8), intent(in), dimension(ni, nj-1) :: vt
  real(8), intent(inout), dimension(ni, nj) :: p
  real(8), dimension(ni,nj) :: p_new, p_old, q

  p_old = p
  p_new = p

  do i=2, ni-1
     do j=2, nj-1
        q(i,j) = (ut(i+1,j) - ut(i,j) + vt(i,j) - vt(i,j-1))/dt*dx
     enddo
  enddo

  do 
     do i= 2, ni-1
        do j= 2, nj-1
           p_new(i,j) = 0.25*beta*(p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - q(i,j)) + (1.0-beta)*p(i,j)
        enddo
     enddo
!     print *, maxval(abs(p_new-p))
     if (maxval(abs(p_new-p)) < tol) then
        p = p_new
        exit
     else 
        p = p_new
        ! do i=1, ni
        !    p(i,1) = p_new(i,2)
        !    p(i,nj) = p_new(i,nj-1)
        ! enddo
        ! do j=1, nj
        !    p(1,j) = p_new(2,j)
        !    p(ni,j) = p_new(ni-1,j)
        ! enddo
     endif
  enddo
end subroutine poisson

subroutine set_boundary(u, v, ut, vt, uw, ni, nj)
  implicit none
  integer ni, nj, i, j
  real(8), intent(inout), dimension(ni+1, nj) :: u, ut
  real(8), intent(inout), dimension(ni, nj-1) :: v, vt
  real(8), intent(in), dimension(4) :: uw
  
  do i=1, ni+1
     do j=1, nj
        if (j.eq.1) then
           u(i,j) = 2*uw(1)-u(i,2)
!           ut(i,j) = 2*uw(1)-ut(i,2)
        else if(j.eq.nj) then
           u(i,j) = 2*uw(2)-u(i,nj-1)
 !          ut(i,j) = 2*uw(2)-ut(i,nj-1)
        endif

        if (i.eq.1) then
           u(i,j) = -u(2,j)
  !         ut(i,j) = -ut(2,j)
        else if(i.eq.ni+1) then
           u(i,j) = -u(ni,j)
   !        ut(i,j) = -ut(ni,j)
        endif
     enddo
  enddo

  do i=1, ni
     do j=1, nj-1
        if (j.eq.1) then
           v(i,j) = -v(i,2)
    !       vt(i,j) = -vt(i,2)
        else if(j.eq.nj-1) then
           v(i,j) = -v(i,nj-2)
     !      vt(i,j) = -vt(i,nj-1)
        endif

        if (i.eq.1) then
           v(i,j) = 2*uw(3)-v(2,j)
      !     vt(i,j) = 2*uw(3)-vt(2,j)
        else if(i.eq.ni) then
           v(i,j) = 2*uw(4)-v(ni-1,j)
       !    vt(i,j) = 2*uw(4)-vt(ni-1,j)
        endif
     enddo
  enddo
     
end subroutine set_boundary

subroutine set_uv_t(u, v, ut, vt, ni, nj, Re, dx, dt)
  implicit none
  integer ni, nj, i, j
  real(8) :: Re, dx, dy, dt
  real(8), intent(in), dimension(ni+1, nj) :: u
  real(8), intent(in), dimension(ni, nj-1) :: v
  real(8), intent(inout), dimension(ni+1, nj) :: ut
  real(8), intent(inout), dimension(ni, nj-1) :: vt
  
  real(8) :: a_1, a_2, a_3, a_4, d_1, d_2, A, D
  dy = dx
  do i=2, ni+1-1
     do j=2, nj-1
        a_1 = (0.5*(u(i+1,j) + u(i,j)))**2.0
        a_2 = (0.5*(u(i,j) + u(i-1,j)))**2.0
        a_3 = 0.25*(u(i,j) + u(i,j+1))*(v(i-1,j) + v(i,j))
        a_4 = 0.25*(u(i,j) + u(i,j-1))*(v(i-1,j-1) + v(i,j-1))
        
        d_1 = (u(i+1,j) - 2.0*u(i,j) + u(i-1,j))/dx/dx
        d_2 = (u(i,j+1) - 2.0*u(i,j) + u(i,j-1))/dy/dy
        
        A = (a_1 - a_2)/dx + (a_3 - a_4)/dy
        D = 1/Re*(d_1 + d_2)
        ut(i,j) = u(i,j) + dt*(-A+D)
     enddo
  enddo
  
  do i=2, ni-1
     do j=2,nj-1-1
        a_1 = (0.5*(v(i,j+1) + v(i,j)))**2.0
        a_2 = (0.5*(v(i,j) + v(i,j-1)))**2.0
        a_3 = 0.25*(u(i+1,j) + u(i+1,j+1))*(v(i,j) + v(i+1,j))
        a_4 = 0.25*(u(i,j) + u(i,j+1))*(v(i-1,j) + v(i,j))
        
        d_1 = (v(i+1,j) - 2.0*v(i,j) + v(i-1,j))/dx/dx
        d_2 = (v(i,j+1) - 2.0*v(i,j) + v(i,j-1))/dy/dy
        
        A = (a_1 - a_2)/dy + (a_3 - a_4)/dx
        D = 1/Re*(d_1 + d_2)
        vt(i,j) = v(i,j) +dt*(-A+D)
     enddo
  enddo
end subroutine set_uv_t

subroutine update_uv(u, v, ut, vt, p, dt, dx, ni, nj)
  implicit none
  integer ni, nj, i, j
  real(8) :: dx, dy, dt
  real(8), intent(inout), dimension(ni+1, nj) :: u
  real(8), intent(inout), dimension(ni, nj-1) :: v
  real(8), intent(in), dimension(ni+1, nj) :: ut
  real(8), intent(in), dimension(ni, nj-1) :: vt
  real(8), intent(in), dimension(ni, nj) :: p

  dy = dx
  do i=2, ni+1-1
     do j=2, nj-1
        u(i,j) = ut(i,j) - dt/dx*(p(i,j)-p(i-1,j))
     enddo
  enddo
  
  do i=2, ni-1
     do j=2,nj-1-1
        v(i,j) = vt(i,j) - dt/dy*(p(i,j+1)-p(i,j))
     enddo
  enddo
  
end subroutine update_uv

subroutine calc_vorticity(u, v, uu, vv, w, dx, dy, ni, nj, n)
  implicit none
  integer ni, nj, n, i, j, ii, jj
  real(8) :: uy1, uy2, vx1, vx2
  real(8), intent(in) :: dx, dy
  real(8), intent(inout), dimension(n,n) :: uu, vv, w
  real(8), intent(in), dimension(ni+1,nj) :: u
  real(8), intent(in), dimension(ni,nj-1) :: v
  do ii=1, n
     do jj=1, n
        i = ii + 1
        j = jj + 1
        uu(ii,jj) = 0.5*(u(i+1,j) + u(i,j))
        vv(ii,jj) = 0.5*(v(i,j) + v(i,j-1))
        w(ii,jj) = (vv(ii+1,jj) - vv(ii,jj))/dx - (uu(ii,jj+1) - uu(ii,jj))/dy
     enddo
  enddo

end subroutine calc_vorticity
   
