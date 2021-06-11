library(DependenciesGraphs)
library(GREAT)

# Prepare data
# dep <- envirDependencies("package:GREAT")
# # dep <- funDependencies("package:GREAT", "shuffle_ro18_timepoints")
#
# # visualization
# plot(dep)

plot(funDependencies("package:GREAT", "get_best_stretch_and_shift"))
