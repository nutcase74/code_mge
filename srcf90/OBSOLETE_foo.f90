!C ---i = 1
      ss = b%si(2,j,k,:)
      fgx1(1,j,k,:) = 0.5d0*(w(1,j,k,1:3)+w(2,j,k,1:3))*ss(1)
      fgy1(1,j,k,:) = 0.5d0*(w(1,j,k,1:3)+w(2,j,k,1:3))*ss(2)
      fgz1(1,j,k,:) = 0.5d0*(w(1,j,k,1:3)+w(2,j,k,1:3))*ss(3)
      fgt1(1,j,k,:) = 0.5d0*(t(1,j,k    )+t(2,j,k    ))*ss

!C ---i = ib
      ss = b%si(ie,j,k,:)
      fgx1(ie,j,k,:) = 0.5d0*(w(il,j,k,1:3)+w(ie,j,k,1:3))*ss(1)
      fgy1(ie,j,k,:) = 0.5d0*(w(il,j,k,1:3)+w(ie,j,k,1:3))*ss(2)
      fgz1(ie,j,k,:) = 0.5d0*(w(il,j,k,1:3)+w(ie,j,k,1:3))*ss(3)
      fgt1(ie,j,k,:) = 0.5d0*(t(il,j,k    )+t(ie,j,k    ))*ss
      
      ss  = 0.5d0*b%sj(il,j,k,:)
      uav = 0.25d0*(w(ie,j  ,k,1:3)+w(il,j  ,k,1:3) &
                   +w(il,j-1,k,1:3)+w(ie,j-1,k,1:3))
      tav = 0.25d0*(t(ie,j  ,k)+t(il,j  ,k) &
                   +t(il,j-1,k)+t(ie,j-1,k))
      fgx2(ie,j,k,:) = uav*ss(1)
      fgy2(ie,j,k,:) = uav*ss(2)
      fgz2(ie,j,k,:) = uav*ss(3)
      fgt2(ie,j,k,:) = tav*ss

      ss  = 0.5d0*b%sk(il,j,k,:)
      uav = 0.25d0*(w(ie,j,k  ,1:3)+w(il,j,k  ,1:3) &
                   +w(ie,j,k-1,1:3)+w(il,j,k-1,1:3))
      tav = 0.25d0*(t(ie,j,k  )+t(il,j,k  ) &
                   +t(ie,j,k-1)+t(il,j,k-1))
      fgx3(ie,j,k,:) = uav*ss(1)
      fgy3(ie,j,k,:) = uav*ss(2)
      fgz3(ie,j,k,:) = uav*ss(3)
      fgt3(ie,j,k,:) = tav*ss




!C ---j = 1
      ss = b%sj(i,2,k,:)
      fgx2(i,1,k,:) = 0.5d0*(w(i,1,k,1:3)+w(i,2,k,1:3))*ss(1)
      fgy2(i,1,k,:) = 0.5d0*(w(i,1,k,1:3)+w(i,2,k,1:3))*ss(2)
      fgz2(i,1,k,:) = 0.5d0*(w(i,1,k,1:3)+w(i,2,k,1:3))*ss(3)
      fgt2(i,1,k,:) = 0.5d0*(t(i,1,k)+t(i,2,k))*ss

!C ---j = je
      ss  = 0.5d0*b%si(i,jl,k,:)
      uav = 0.25d0*(w(i,je,k,1:3)+w(i-1,je,k,1:3) &
                   +w(i,jl,k,1:3)+w(i-1,jl,k,1:3))
      tav = 0.25d0*(t(i,je,k)+t(i-1,je,k) &
                   +t(i,jl,k)+t(i-1,jl,k))
      fgx1(i,je,k,:) = uav*ss(1)
      fgy1(i,je,k,:) = uav*ss(2)
      fgz1(i,je,k,:) = uav*ss(3)
      fgt1(i,je,k,:) = tav*ss
      
      ss = b%sj(i,je,k,:) 
      fgx2(i,je,k,:) = 0.5d0*(w(i,je,k,1:3)+w(i,jl,k,1:3))*ss(1)
      fgy2(i,je,k,:) = 0.5d0*(w(i,je,k,1:3)+w(i,jl,k,1:3))*ss(2)
      fgz2(i,je,k,:) = 0.5d0*(w(i,je,k,1:3)+w(i,jl,k,1:3))*ss(3)
      fgt2(i,je,k,:) = 0.5d0*(t(i,je,k)+t(i,jl,k))*ss
      
      ss = 0.5d0*b%sk(i,jl,k,:)
      uav = 0.25d0*(w(i,je,k  ,1:3)+w(i,jl,k  ,1:3) &
                   +w(i,je,k-1,1:3)+w(i,jl,k-1,1:3))
      tav = 0.25d0*(t(i,je,k  )+t(i,jl,k  ) &
                   +t(i,je,k-1)+t(i,jl,k-1))
      fgx3(i,je,k,:) = uav*ss(1)
      fgy3(i,je,k,:) = uav*ss(2)
      fgz3(i,je,k,:) = uav*ss(3)
      fgt3(i,je,k,:) = tav*ss



!C ---k = 1
      ss = b%sk(i,j,2,:)
      fgx3(i,j,1,:) = 0.5d0*(w(i,j,1,1:3)+w(i,j,2,1:3))*ss(1)
      fgy3(i,j,1,:) = 0.5d0*(w(i,j,1,1:3)+w(i,j,2,1:3))*ss(2)
      fgz3(i,j,1,:) = 0.5d0*(w(i,j,1,1:3)+w(i,j,2,1:3))*ss(3)
      fgt3(i,j,1,:) = 0.5d0*(t(i,j,1)+t(i,j,2))*ss

!C ---k = ke
      ss  = 0.5d0*b%si(i,j,kl,:)
      uav = 0.25d0*(w(i  ,j,ke,1:3)+w(i-1,j,ke,1:3) &
                   +w(i-1,j,kl,1:3)+w(i  ,j,kl,1:3))
      tav = 0.25d0*(t(i  ,j,ke)+t(i-1,j,ke) &
                   +t(i-1,j,kl)+t(i  ,j,kl))
      fgx1(i,j,ke,:) = uav*ss(1)
      fgy1(i,j,ke,:) = uav*ss(2)
      fgz1(i,j,ke,:) = uav*ss(3)
      fgt1(i,j,ke,:) = tav*ss

      ss  = 0.5d0*b%sj(i,j,kl,:)
      uav = 0.25d0*(w(i,j  ,ke,1:3)+w(i,j  ,kl,1:3) &
                   +w(i,j-1,ke,1:3)+w(i,j-1,kl,1:3))
      tav = 0.25d0*(t(i,j  ,ke)+t(i,j  ,kl) &
                   +t(i,j-1,ke)+t(i,j-1,kl))
      fgx2(i,j,ke,:) = uav*ss(1)
      fgy2(i,j,ke,:) = uav*ss(2)
      fgz2(i,j,ke,:) = uav*ss(3)
      fgt2(i,j,ke,:) = tav*ss

      ss = b%sk(i,j,ke,:)
      fgx3(i,j,ke,:) = 0.5d0*(w(i,j,ke,1:3)+w(i,j,kl,1:3))*ss(1)
      fgy3(i,j,ke,:) = 0.5d0*(w(i,j,ke,1:3)+w(i,j,kl,1:3))*ss(2)
      fgz3(i,j,ke,:) = 0.5d0*(w(i,j,ke,1:3)+w(i,j,kl,1:3))*ss(3)
      fgt3(i,j,ke,:) = 0.5d0*(t(i,j,ke)+t(i,j,kl))*ss
