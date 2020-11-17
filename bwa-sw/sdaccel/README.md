Get datafiles
-------------
%cd testdata  
%./get-data.sh  
%cd ..  

Run pure C simulation
---------------------
%make clean  
%make run_csim  

-----------------------------------------------------
For targetting non-F1 FPGA
-----------------------------------------------------
Set up environment
------------------
%module load sdx/16.3  

Run software emulation
----------------------
%make clean  
%make run_sw_emu  

Run hardware emulation
----------------------
%make clean  
%make run_hw_emu  

Generate kernel binary
----------------------
%make clean  
%make bit  

Run tests on hardware
--------------------
%make clean  
%make run  

-----------------------------------------------------
For targetting F1 FPGA 
-----------------------------------------------------
Set up environment
------------------
%module load sdx/17.1  

Run software emulation
----------------------
%make clean  
%make XDEVICE=F1 run_sw_emu  

Run hardware emulation
----------------------
%make clean  
%make XDEVICE=F1 run_hw_emu  

Generate kernel binary
----------------------
%make clean  
%make XDEVICE=F1 bit  

Run tests on hardware
--------------------
%make clean  
%make XDEVICE=F1 run  

