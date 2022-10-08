# insulation score of MC and FL at Sox9 locus
cooler coarsen -k 5 -n 5 -o MC_HiC_combined_50kb.cool /data05/CQM/CQM_20200401_HiC_combined/E145_MC_combined_10kb.cool
cooler coarsen -k 5 -n 5 -o FL_HiC_combined_50kb.cool /data05/CQM/CQM_20200401_HiC_combined/E145_Forelimb_combined_10kb.cool
cooler balance MC_HiC_combined_50kb.cool
cooler balance FL_HiC_combined_50kb.cool
cooltools insulation  MC_HiC_combined_50kb.cool --window-pixels 10 -o  MC_HiC_combined_50kb_insulation --bigwig 
cooltools insulation  FL_HiC_combined_50kb.cool  --window-pixels 10 -o FL_HiC_combined_50kb_insulation --bigwig 

