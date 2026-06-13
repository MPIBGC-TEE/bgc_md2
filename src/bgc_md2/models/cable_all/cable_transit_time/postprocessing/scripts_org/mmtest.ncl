begin
DirPathStr="/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4/"
ntime=365
npatch=10
nland=5656
npool=9
B=new((/ntime,npool,npool,npatch,nland/),float)
C=new((/npool,npool,npatch,nland/),float)
C=0
Ct=C
do day=0,ntime-1
	Ct:=array_append_record(Ct,Ct,0)
end do
year = 1901
	fin = addfile(DirPathStr+"/out_ncar_"+year+"_ndep.nc","r")
	iveg_org=fin->iveg
	do day=0,ntime-1
		do ipool=0,npool-1
			do jpool=0,npool-1
				B(day,ipool,jpool,:,:)=onedtond(\
		      			 where(\
						ismissing(iveg_org),\
						B@_FillValue,\
						0\
					 ),\
					 (/npatch,nland/)\
				) 
			end do
		end do
	end do
	;B(:,3,1,:,:) = fin->fromRoottoL(:,0,:,:)*fin->kplant(:,2,:,:)
	;f. Metabolic turnover
	x= Ct(:3,:,:)*fin->xktemp(:,:,:)
print(B)
end
