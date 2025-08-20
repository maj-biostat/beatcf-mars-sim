

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

report rep:
  quarto render reports/{{rep}} --to pdf



