## Analysis of copy number structural variants by pair-end short read simulation

### Introduction
Structural variations (SV), in particular copy number variations, play an important role in serious health conditions, such as various cancers. Presence and severity of copy number variations are in direct correlation with the condition and progress of cancerous cells. State-of-the-art SV callers, bioinformatics tools that detect SVs from sequencing data, do not effectively detect and characterize copy number variations. This study was designed to simulation RNA-seq reads that contain pre-defined copy number variations to assess novel detection algorithm design from the sequencing data.

### Scripts
`segmentCN.R`: Segmentation library for copy number assessment of a given region.
`CN_profiles.Rmd`: Provide show cases of interesting copy number variations in the simulated data.
`Simulation_functions.R`: RNA-seq data simulation library
`thesis_results.Rmd`: Contains a subset of my master's thesis results.
