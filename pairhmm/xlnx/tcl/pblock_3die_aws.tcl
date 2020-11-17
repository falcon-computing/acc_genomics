set SLR0_pblock pblock_CL_bot 
set SLR1_pblock pblock_CL_mid
set SLR2_pblock pblock_CL_top 


add_cells_to_pblock $SLR0_pblock [get_cells -hier *pmm_core_top0*] -clear_locs
add_cells_to_pblock $SLR1_pblock [get_cells -hier *pmm_core_top1*] -clear_locs
add_cells_to_pblock $SLR2_pblock [get_cells -hier *pmm_core_top2*] -clear_locs

