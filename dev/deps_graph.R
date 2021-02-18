library(DependenciesGraphs)
library(genereg)

# Prepare data
# dep <- envirDependencies("package:genereg")
# # dep <- funDependencies("package:genereg", "shuffle_ro18_timepoints")
#
# # visualization
# plot(dep)

plot(funDependencies("package:genereg", "get_best_stretch_and_shift"))
