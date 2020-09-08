
:: This batch file checks input _tri.raw,
@echo on
:: OUTPUT
:: OUTPUT
@echo Run segmentation
Segmentation.exe -i rod_tri.k -o rod_initial_write.k -m rod_manual.txt -l 0.1
@echo Done!
@echo -------------------------------------------------------------------
@pause
:: OUTPUT
@echo Create polycube structure
PolyCube.exe -i rod_initial_read.k -o rod_polycube_structure.k -c 1
@echo Done!
@echo -------------------------------------------------------------------
@pause
:: OUTPUT
@echo Generate hex mesh
ParametricMapping.exe -i rod_initial_read.k -p rod_polycube_structure_hex.k -o rod_hex.vtk -s 2
@echo Done!
@echo -------------------------------------------------------------------
@pause
@echo Fix interior points 
Quality.exe -I rod_hex.vtk   -m 0 -n 1 
@echo Done!
@echo -------------------------------------------------------------------
@pause

@echo Pillowing
Quality.exe -I rod_hex_lap.vtk   --method 1 --number 1
@echo Done!
@echo If needed, prepare sharp feature in sharp.txt before moving to next step
@echo -------------------------------------------------------------------
@pause

@echo Smooth  points
Quality.exe -I ./rod_hex_lap_pillow.vtk --method 2 --parameter 0.01 --number 1000 --sharp 2 
@echo Done!
@echo -------------------------------------------------------------------
@pause

@echo Optimize element with minimum Jacobian
Quality.exe -I rod_hex_lap_pillow_smooth.vtk   --method 3 --parameter 0.01 --number 150
@echo Done!
@echo -------------------------------------------------------------------
@pause

@echo Spline construction of Rod model with no refinement and output LS-DYNA simulation file
Hex2Spline.exe -I rod_hex.vtk -S -s 2 
@echo Done!
@echo -------------------------------------------------------------------
@pause

@echo Spline construction of Rod model with one global refinement and output LS-DYNA simulation file
Hex2Spline.exe -I rod_hex.vtk -S -s 2 -g 1
@echo Done!
@echo -------------------------------------------------------------------
@pause

@echo Spline construction of Rod model with level-2 local refinement and output LS-DYNA simulation file
Hex2Spline.exe -I rod_hex.vtk -S -s 2 -l 
@echo Done!
@echo -------------------------------------------------------------------
@pause

 
