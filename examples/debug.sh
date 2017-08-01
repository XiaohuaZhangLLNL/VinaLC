#!/bin/bash

srun -N4 -n4 -c12 -ppdebug vinaBMPI --recList recList.txt --ligList ligList.txt --geoList geoList.txt 
#srun -N4 -n4 -c12 -ppdebug ./vina --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/NetBeansProjects/VinaLC/apps/vinaBMPI --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 --num_modes 20
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/bin/vinaBMPI --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 --num_modes 20
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/NetBeansProjects/VinaLC/apps/vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 --num_modes 20
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/NetBeansProjects/VinaLC/apps/vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 
