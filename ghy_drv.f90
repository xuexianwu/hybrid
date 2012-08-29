!======================================================================!
subroutine ghinij
!----------------------------------------------------------------------!
! For local position:
! calculate number of soil layers (n);
!----------------------------------------------------------------------!
use control_parameters
use constants
use soil_variables
use veg_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: i ! Soil texture index                                    (n)
!----------------------------------------------------------------------!
! Get local soil layer thicknesses (m).
!----------------------------------------------------------------------!
dz (1:nsoil_layers_max) = dz_ij (i0,j0,1:nsoil_layers_max)
!----------------------------------------------------------------------!
! Get local soil texture fractions.
!----------------------------------------------------------------------!
q (1:ntextures,1:nsoil_layers_max) = &
& q_ij (i0,j0,1:ntextures,1:nsoil_layers_max)
!----------------------------------------------------------------------!
! Get local soil texture fractions (for conductivities).
!----------------------------------------------------------------------!
qk (1:ntextures,1:nsoil_layers_max) = &
& qk_ij (i0,j0,1:ntextures,1:nsoil_layers_max)
!----------------------------------------------------------------------!
! Local slope (units?).
!----------------------------------------------------------------------!
sl = sl_ij (i0,j0)
!----------------------------------------------------------------------!
! Calculate number of soil layers.
!----------------------------------------------------------------------!
nsoil_layers = 0
do nsoil_layer = 1, nsoil_layers_max
  if (dz (nsoil_layer) <= 0.0) exit
  nsoil_layers = nsoil_layer
enddo
if (nsoil_layers <= 0) then
  write (*,*) 'no soil layers',i0,j0
  write (99,*) 'no soil layers',i0,j0
  stop
endif
!----------------------------------------------------------------------!
! Calculate boundaries, based on thicknesses.
!----------------------------------------------------------------------!
zb (1) = 0.0d0
do nsoil_layer = 1, nsoil_layers
  zb (nsoil_layer+1) = zb (nsoil_layer) - dz (nsoil_layer)
enddo
!----------------------------------------------------------------------!
! Calculate layer centres, based on boundaries.
!----------------------------------------------------------------------!
do nsoil_layer = 1, nsoil_layers
  zc (nsoil_layer) = 0.5d0 * (zb (nsoil_layer) + zb (nsoil_layer+1))
enddo
fb = afb (i0,j0)
fv = 1.0d0 - fb
!----------------------------------------------------------------------!
! Calculate minimum and saturated thetas (i.e. volumetric moisture
! content) for each layer (m3/m3), and saturated water content (m).
!----------------------------------------------------------------------!
do ibv = 1, 2
  do nsoil_layer = 1, nsoil_layers
    thets (nsoil_layer,ibv) = 0.0
    thetm (nsoil_layer,ibv) = 0.0
    do i = 1, ntextures - 1
      thets (nsoil_layer,ibv) = thets (nsoil_layer,ibv) + &
      &                         q (i,nsoil_layer) * thm (0,i)
      thetm (nsoil_layer,ibv) = thetm (nsoil_layer,ibv) + &
      &                         q (i,nsoil_layer) * thm (nth,i)
    enddo
    ws (nsoil_layer,ibv) = thets (nsoil_layer,ibv) * dz (nsoil_layer)
  enddo
enddo
!----------------------------------------------------------------------!
! Calculate specific heat capacities of layers (J/(m^2 deg C))
!----------------------------------------------------------------------!
do ibv = 1, 2
  do nsoil_layer = 1, nsoil_layers
    shc (nsoil_layer,ibv) = 0.0d0
    do i = 1, ntextures ! Loop over textures.
      shc (nsoil_layer,ibv) = shc (nsoil_layer,ibv) + &
      &                       q (i,nsoil_layer) * shcap (i)
    enddo
    shc (nsoil_layer,ibv) = (1.0d0 - thets (nsoil_layer,ibv)) * &
    &                       shc (nsoil_layer,ibv) * dz (nsoil_layer)
  enddo
enddo
!----------------------------------------------------------------------!
end subroutine ghinij
!======================================================================!
