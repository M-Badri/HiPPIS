module machine

integer, parameter :: r8_fix    = SELECTED_REAL_KIND(12)   ! real r8
integer, parameter :: r4_fix    = SELECTED_REAL_KIND(6,30) ! real r4
integer, parameter :: real_kind = r8_fix
integer, parameter :: kind_phys = real_kind

end module
