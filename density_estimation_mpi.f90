module fun
contains

subroutine kernel_density_estimate(m,y_s,n,z_s,points,bandwidth,map)
        implicit none
        real :: map
        real :: points(:,:)
        integer :: s
        real :: bandwidth
        real,allocatable :: distance(:)
        real,allocatable :: kernel_values(:)
        integer :: m,n
        integer :: k
        real,parameter :: pi=3.14159
        real :: y_s,z_s

        s=size(points(1,:))
        allocate(distance(s))
        allocate(kernel_values(s))

        do k=1,s
                distance(k)=(points(1,k)-m*y_s)**2+(points(2,k)-n*z_s)**2
                kernel_values(k)=exp(-0.5*distance(k)/(bandwidth**2))
        end do

        map=sum(kernel_values)/((bandwidth**2))/s
        deallocate(distance)
        deallocate(kernel_values)
end subroutine kernel_density_estimate
end module fun

program density_estimation
        use fun
        !$ use omp_lib
        implicit none
        type :: atom_3d
                integer :: atom_type
                real :: x(3)
        end type atom_3d
        real,allocatable :: x_b(:,:),a(:),b(:),c(:)
        integer,allocatable :: atom_number(:,:)
        real,allocatable :: density_map(:,:,:),density_map_average(:,:)
        integer :: i,j,k,l,step,alive,n,ierr,num_processors,my_processor_id
        type(atom_3d),allocatable :: R_3d(:,:)
        real,allocatable :: R_2d(:,:,:)
        real,allocatable :: y_step(:),z_step(:),y(:)
        real :: bandwidth,density_max,density_min,y_step_a,z_step_a
        integer,parameter :: n_points=1000
        integer,parameter :: cores=10
        character(len=1024) :: step_i,read_file
!-----------------------------------------------------------------
write(*,*) "please write down 1. read_file  2. total step"
read(*,*) read_file,step_i
read(step_i,*) step
open(unit=100,file=read_file,status='old',action='read')
open(unit=200,file='O/density-x.dat',status='replace')
open(unit=300,file='O/density-y.dat',status='replace')
open(unit=400,file='O/density-z.dat',status='replace')
open(unit=600,file='O/projection.dat',status='replace')
!----------------------------------------------------------------
allocate(atom_number(5,step))
!----------------------------------------------------------------
read(100,*)
read(100,*)
read(100,*)
read(100,*) atom_number(1,1) 
rewind(100)
!----------------------------------------------------------------
allocate(R_3d(atom_number(1,1),step))
allocate(density_map(n_points,n_points,step))
allocate(density_map_average(n_points,n_points))
allocate(x_b(6,step))
allocate(a(step))
allocate(b(step))
allocate(c(step))
allocate(y_step(step))
allocate(z_step(step))
allocate(y(n_points))
!----------------------------------------------------------------
n=1
atom_number(2:5,:)=0
do while (.true.)
        read(100,*,iostat=alive)
        if (alive/=0) then
                goto 55
        end if
        read(100,*)
        read(100,*)
        read(100,*) atom_number(1,n)
        read(100,*)
        read(100,*) x_b(1,n),x_b(2,n)
        read(100,*) x_b(3,n),x_b(4,n)
        read(100,*) x_b(5,n),x_b(6,n)
        a(n)=x_b(2,n)-x_b(1,n)
        b(n)=x_b(4,n)-x_b(3,n)
        c(n)=x_b(6,n)-x_b(5,n)
        y_step(n)=b(n)/n_points
        z_step(n)=c(n)/n_points
        read(100,*)

        do i=1,atom_number(1,n)
                read(100,*) l,R_3d(l,n)%atom_type,R_3d(l,n)%x(1),R_3d(l,n)%x(2),R_3d(l,n)%x(3)
                R_3d(l,n)%x(1)=R_3d(l,n)%x(1)+abs(x_b(1,n))
                R_3d(l,n)%x(2)=R_3d(l,n)%x(2)+abs(x_b(3,n))
                R_3d(l,n)%x(3)=R_3d(l,n)%x(3)+abs(x_b(5,n))
                if (R_3d(l,n)%atom_type == 1) then
                        atom_number(2,n)=atom_number(2,n)+1
                else if (R_3d(l,n)%atom_type == 2) then
                        atom_number(3,n)=atom_number(3,n)+1
                else if (R_3d(l,n)%atom_type == 3) then
                        atom_number(4,n)=atom_number(4,n)+1
                else if (R_3d(l,n)%atom_type == 4) then
                        atom_number(5,n)=atom_number(5,n)+1
                end if
        end do
        n=n+1
end do
55 close(100)
deallocate(x_b)
deallocate(a)
deallocate(b)
deallocate(c)
allocate(R_2d(2,atom_number(5,1),step)) ! if you changed the element, here is a must
write(*,*) "reading finished!"
!--------------------------------------------------------------------------------------------
do i=1,step
        write(600,*) "step:",(i-1)*100
        l=1
loop:   do j=1,atom_number(1,i)
                if (R_3d(j,i)%atom_type == 4) then
                        R_2d(1,l,i)=R_3d(j,i)%x(2) ! here we make a projection.
                        R_2d(2,l,i)=R_3d(j,i)%x(3)
                        write(600,*) R_2d(1,l,i),R_2d(2,l,i)
                        l=l+1
                else
                        cycle loop
                end if
        end do loop
end do
deallocate(R_3d)
write(*,*) "projection finished"
close(600)
!--------------------------------------------------------------------------------------------
bandwidth=0.1
density_map_average(:,:)=0.0
density_map(:,:,:)=0.0
!$omp parallel num_threads(cores) private(n)
!$ n=omp_get_thread_num()
do j=(n_points/cores)*n+1,(n_points/cores)*(n+1)
        do k=1,n_points
                do i=1,step
                        call kernel_density_estimate(j,y_step(i),k,z_step(i),R_2d(:,:,i),bandwidth,density_map(j,k,i))
                        density_map_average(j,k)=density_map_average(j,k)+density_map(j,k,i)
                end do
                density_map_average(j,k)=density_map_average(j,k)/(1.0*step)
        end do
        write(*,*) 'points:',j,k,"finished!"
end do
!$omp end parallel
!--------------------------------------------------------------------------------------------------
density_max=maxval(density_map_average)
density_min=minval(density_map_average)
do i=1,n_points
        do j=1,n_points
                density_map_average(i,j)=(density_map_average(i,j)-density_min)/(density_max-density_min)
        end do
end do

y_step_a=0.0
z_step_a=0.0
do i=1,step
        y_step_a=y_step(i)+y_step_a
        z_step_a=z_step(i)+z_step_a
end do
deallocate(y_step)
deallocate(z_step)
y_step_a=y_step_a/(1.0*step)
z_step_a=z_step_a/(1.0*step)
do i=1,n_points
        write(200,*) (j*y_step_a,j=1,n_points)
        y(:)=i*z_step_a
        write(300,*) (y(k),k=1,n_points)
        write(400,*) (density_map_average(j,i),j=1,n_points)
end do
deallocate(R_2d)
deallocate(density_map)
stop
end program
