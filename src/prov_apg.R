library(provR)
library(Rclean)
prov.capture("apg.R")
tstamp <- format(Sys.time(), "%H_%d_%m_%Y")
file.out <- paste0("../data/storage/apg/apg_R_", tstamp, ".json")
write.code(prov.json(), file.out)


