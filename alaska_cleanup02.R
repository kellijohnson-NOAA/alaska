
# unload the cpp file
dyn.unload(dynlib(file.path(dir.data, my.tmb)))
gc()
try(dyn.unload(dynlib(file.path(dir.data, my.tmb))), silent = TRUE)
