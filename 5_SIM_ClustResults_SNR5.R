#
# Simulation study - clustering performance under low noise
# as reported in Section 3.1 of the main manuscript
#
load("SIM_clusteringResults_SNR5_Inter")
# 
#
######################################################
# Low noise (STN5), eigendimension 1
######################################################
#
# 1) # of simulation correctly classified: ARI vs VI
#
table(round(ARI1,3))
#0.907| 0.951 |    1 
#   3 |    4  |   93
#
table(round(VI_classif,3))
#0.951|     1 
#    3|    97 
which(VI_classif<1)
#
## uncertainty with Cluster balls
#
# Horizontal, upper vertical, lower vertical
View(uncert_VI_h) # 5/97: groupB-cluster->subject-specific cluster
View(uncert_VI_u) # none ( minimal clustering most distant
                  # to point estimate coincides with point estimates)
View(uncert_VI_l) # 5/97: groupB-cluster->subject-specific cluster (same simulations)
#
######################################################
# Low noise (STN5), eigendimension 2
######################################################
#
# 1) # of simulation correctly classified: ARI vs VI
#
table(round(ARI2,3))
#    1 
#   100
#
table(round(VI2_classif,3))
#    1 
#   100 
#
## uncertainty with Cluster balls
#
# Horizontal, upper vertical, lower vertical
View(uncert_VI2_h) # 1/97: groupB-cluster->subject-specific cluster
View(uncert_VI2_u) # none (minimal clustering most distant
                   # to point estimate coincides with point estimates)
View(uncert_VI2_l) # 1/97: groupB-cluster->subject-specific cluster (same simulation)
#
## clustering in cluster-specific clusters (u=1,2,39,40)
#
table(round(VI_s1_res,2))
#0.92|    1 
#  14|   86 
#
table(round(VI_s2_res,2))
#0.77|  0.84|  0.92|     1 
#  34|    16|    17|    33 
#
table(round(VI_s39_res,2))
#0.92|    1 
#   5|   95 
#
table(round(VI_s40_res,2))
#0.77| 0.84 
#   3|   97
#
# 63% of all individual deviations classified
# with 0-1 errors/50 (ARI>0.9) 
# 37% of all individual deviations classified
# with 2-3 errors/50 (ARI>0.7)
#
