#!/bin/bash

ncremap -P mpas -i subgrid_on_subgridout_test_3.2000-01-02.nc -m /home/zwolff/Models/MPAS_SeaIce/Output_test/map_QU120_seaice_to_r05_202101.nc   -o subgrid_on_test_3_2000-01-02.nc 

ncremap -P mpas -i subgrid_on_subgridout_test_3.2000-01-03.nc -m /home/zwolff/Models/MPAS_SeaIce/Output_test/map_QU120_seaice_to_r05_202101.nc   -o subgrid_on_test_3_2000-01-03.nc 

ncremap -P mpas -i subgrid_on_subgridout_test_3.2000-01-04.nc -m /home/zwolff/Models/MPAS_SeaIce/Output_test/map_QU120_seaice_to_r05_202101.nc   -o subgrid_on_test_3_2000-01-04.nc 

ncremap -P mpas -i subgrid_off_subgridout_test_3.2000-01-02.nc -m /home/zwolff/Models/MPAS_SeaIce/Output_test/map_QU120_seaice_to_r05_202101.nc   -o subgrid_off_test_3_2000-01-02.nc 

ncremap -P mpas -i subgrid_off_subgridout_test_3.2000-01-03.nc -m /home/zwolff/Models/MPAS_SeaIce/Output_test/map_QU120_seaice_to_r05_202101.nc   -o subgrid_off_test_3_2000-01-03.nc 

ncremap -P mpas -i subgrid_off_subgridout_test_3.2000-01-04.nc -m /home/zwolff/Models/MPAS_SeaIce/Output_test/map_QU120_seaice_to_r05_202101.nc   -o subgrid_off_test_3_2000-01-04.nc 
