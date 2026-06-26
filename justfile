

site:
  quarto render --to html

cleantmp:
  find tmp -delete
  mkdir tmp


sim00 cfg:
  Rscript --vanilla ./R/sim00.R run_sim00 {{cfg}} 
  just cleantmp

runsim00:
  just sim00 ../etc/sim00/cfg-sim00-v01.yml
  just sim00 ../etc/sim00/cfg-sim00-v02.yml
  just sim00 ../etc/sim00/cfg-sim00-v03.yml
  just sim00 ../etc/sim00/cfg-sim00-v04.yml
  just sim00 ../etc/sim00/cfg-sim00-v05.yml
  just sim00 ../etc/sim00/cfg-sim00-v06.yml
  just sim00 ../etc/sim00/cfg-sim00-v07.yml

sim01 cfg:
  Rscript --vanilla ./R/sim01.R run_sim01 {{cfg}} 
  just cleantmp

runsim01:
  just sim01 ../etc/sim01/cfg-sim01-v01.yml
  just sim01 ../etc/sim01/cfg-sim01-v02.yml
  just sim01 ../etc/sim01/cfg-sim01-v03.yml
  just sim01 ../etc/sim01/cfg-sim01-v04.yml
  just sim01 ../etc/sim01/cfg-sim01-v05.yml
  just sim01 ../etc/sim01/cfg-sim01-v06.yml
  just sim01 ../etc/sim01/cfg-sim01-v07.yml

sim02 cfg:
  Rscript --vanilla ./R/sim02.R run_sim02 {{cfg}}
  just cleantmp

runsim02:
  just sim02 ../etc/sim02/cfg-sim02-v01.yml
  just sim02 ../etc/sim02/cfg-sim02-v02.yml
  just sim02 ../etc/sim02/cfg-sim02-v03.yml

sim03 cfg:
  Rscript --vanilla ./R/sim03.R run_sim03 {{cfg}}
  just cleantmp

runsim03:
  just sim03 ../etc/sim03/cfg-sim03-v01.yml
  just sim03 ../etc/sim03/cfg-sim03-v02.yml
  just sim03 ../etc/sim03/cfg-sim03-v03.yml
  just sim03 ../etc/sim03/cfg-sim03-v04.yml

sim04 cfg:
  Rscript --vanilla ./R/sim04.R run_sim04 {{cfg}}
  just cleantmp

runsim04:
  just sim04 ../etc/sim04/cfg-sim04-v01.yml
  just sim04 ../etc/sim04/cfg-sim04-v02.yml
  just sim04 ../etc/sim04/cfg-sim04-v03.yml
  just sim04 ../etc/sim04/cfg-sim04-v04.yml
  

sim12 cfg:
  Rscript --vanilla ./R/sim12.R run_sim12 {{cfg}}
  just cleantmp

runsim12:
  just sim12 ../etc/sim12/cfg-sim12-v01.yml
  just sim12 ../etc/sim12/cfg-sim12-v02.yml
  just sim12 ../etc/sim12/cfg-sim12-v03.yml
  

sim14 cfg:
  Rscript --vanilla ./R/sim14.R run_sim14 {{cfg}}
  just cleantmp

runsim14:
  #just sim14 ../etc/sim14/cfg-sim14-v01.yml
  just sim14 ../etc/sim14/cfg-sim14-v02.yml
  #just sim14 ../etc/sim14/cfg-sim14-v03.yml
  
sim15 cfg:
  Rscript --vanilla ./R/sim15.R run_sim15 {{cfg}}
  just cleantmp

runsim15:
  just sim15 ../etc/sim15/cfg-sim15-v01.yml
  #just sim15 ../etc/sim15/cfg-sim15-v02.yml
  
sim16 cfg:
  Rscript --vanilla ./R/sim16.R run_sim16 {{cfg}}
  just cleantmp

runsim16:
  just sim16 ../etc/sim16/cfg-sim16-v01.yml  
  #just sim16 ../etc/sim16/cfg-sim16-v02.yml  
  #just sim16 ../etc/sim16/cfg-sim16-v03.yml  
  #just sim16 ../etc/sim16/cfg-sim16-v04.yml  
  
sim18 cfg:
  Rscript --vanilla ./R/sim18.R run_sim18 {{cfg}}
  just cleantmp

runsim18:
  just sim18 ../etc/sim18/cfg-sim18-v01.yml  
  #just sim18 ../etc/sim18/cfg-sim18-v02.yml  
  
report rep:
  quarto render reports/{{rep}} --to pdf



