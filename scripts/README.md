This folder contains scripts, environments and installation commands to run Coinfinder, Goldfinder, and EvoScope and measure their runtimes.

The commands are given in ´anaysis.sh´, but I did not manage to use mamba within a .sh script - you will have to copy paste the commands into your shell.

I measure the time by creating files start.time and end.time before and after the respective run. Their creation time can be read out using stat ´<file>´. Coinfinder calculates all of the D-values after its execution. So to measure the time in a fair manner, ´stat coin_pairs.tsv´ instead of ´end.time´.

I ran Goldfinder and Coinfinder in this way. Coinfinder took 00:23:30.693435992 (23 min), Goldfinder 00:25:20.195655488 (25 min).