#!/bin/bash

srun -N4 -n4 -c12 -ppdebug vinalc --recList recList.txt --ligList ligList.txt --geoList geoList.txt 
#srun -N4 -n4 -c12 -ppdebug ./vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/NetBeansProjects/VinaLC/apps/vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 --num_modes 20
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/bin/vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 --num_modes 20
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/NetBeansProjects/VinaLC/apps/vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 --num_modes 20
#srun -N4 -n4 -c12 -ppdebug /g/g92/zhang30/medchem/NetBeansProjects/VinaLC/apps/vinalc --recList rigList --ligList ligList.txt --geoList geoList.txt --fleList fleList --randomize 1 
