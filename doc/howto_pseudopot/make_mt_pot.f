        dimension r(1000),v1(1000),v2(1000),v3(1000)
        dimension ps1(1000),ps2(1000),ps3(1000)
        open(unit=1,file='Si_PBE.ps',status='old')
        open(unit=2,file='Si_PBE.wf',status='old')
        open(unit=10,file='Si_l0',status='new')
        open(unit=11,file='Si_l1',status='new')
        open(unit=12,file='Si_l2',status='new')
        read(1,*) npts
        do i=1,npts
         read(1,*) r(i),v1(i),v2(i),v3(i)
        enddo
        read(2,*) npts
        do i=1,npts
         read(2,*) r(i),ps1(i),ps2(i),ps3(i)
        enddo
        write(10,*) npts
        do i=1,npts
         write(10,*) r(i),v1(i),ps1(i)
        enddo
        write(11,*) npts
        do i=1,npts
         write(11,*) r(i),v2(i),ps2(i)
        enddo
        write(12,*) npts
        do i=1,npts
         write(12,*) r(i),v3(i),ps3(i)
        enddo
        stop
        end
        
