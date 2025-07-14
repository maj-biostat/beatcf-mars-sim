

site:
  quarto render --to html

cleantmp:
  find tmp -delete
  mkdir tmp


sim01 cfg:
  Rscript --vanilla ./R/sim01.R run_sim01 {{cfg}} 
  just cleantmp

runsim:
  just sim01 ../etc/sim01/cfg-sim01-v01.yml
  just sim01 ../etc/sim01/cfg-sim01-v02.yml
  just sim01 ../etc/sim01/cfg-sim01-v03.yml
  just sim01 ../etc/sim01/cfg-sim01-v04.yml
  just sim01 ../etc/sim01/cfg-sim01-v05.yml
  just sim01 ../etc/sim01/cfg-sim01-v06.yml
  just sim01 ../etc/sim01/cfg-sim01-v07.yml

report rep:
  quarto render reports/{{rep}} --to pdf



