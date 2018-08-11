module element_Z
! use gauss_rule
! use analytic_part
contains 


subroutine element_Zmn_Full(edge_m,edge_n,value,msh,ker)
    
    use MODULE_FILE
    implicit none
    
    integer edge_m,edge_n
    complex(kind=8) value
    type(mesh)::msh
	type(kernelquant)::ker
	
	value = ker%matZ_glo(msh%new2old(edge_m),msh%new2old(edge_n))
	

end subroutine element_Zmn_Full


    

end module element_Z
    