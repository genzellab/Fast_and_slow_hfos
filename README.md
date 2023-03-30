# Fast and slow cortical high frequency oscillations for cortico-cortical and corticohippocampal network consolidation during NonREM sleep. 

Adrian Aleman-Zapata 1, Richard GM Morris 2, Lisa Genzel *1,2

*corresponding author: lgenzel@donders.ru.nl  :mailbox: 

1, Donders Institute for Brain Cognition and Behavior, Radboud University, Postbus 9010, 6500GL Nijmegen/Netherlands.

2, Centre for Cognitive and Neural Systems, Edinburgh Neuroscience, University of Edinburgh, 1 George Square, Edinburgh EH8 9JZ, UK.

Reference:  [doi.org/10.1101/765149](https://doi.org/10.1101/765149) 

-----------------------------




:warning: Requirements: Makes use of functions from the [YASA](https://github.com/raphaelvallat/yasa) (only used in the older version of the manuscript 2020), [Fieldtrip](https://github.com/fieldtrip/fieldtrip) (Suggest to use 2019 version to avoid plotting issues) and [ADRITOOLS](https://github.com/Aleman-Z/ADRITOOLS) repositories. 

These last two need to be added to the path.

We also include in the subfunctions folder some functions borrowed from [FMA toolbox](https://github.com/michael-zugaro/FMAToolbox/tree/master/Analyses)

_Required Matlab built-in toolboxes:_
•	Image Processing Toolbox
•	Signal Processing Toolbox
•	Statistics and Machine Learning Toolbox
•	Mapping Toolbox
•	Deep Learning Toolbox
•	Symbolic Math Toolbox
•	Bioinformatics Toolbox
•	Computer Vision Toolbox
•	Fixed-Point Designer
•	MATLAB Coder
•	Simulink
•	Parallel Computing Toolbox
•	MATLAB Parallel Server
•	Polyspace Bug Finder

--------------------------------
## Main scripts: :file_folder: 

_**Figure 2B and 2D:** Count of coocurring and single events, as well as slow and fast counts and rates._
  * GL_hfos_counts.m

_Mentioned in Text: Shuffling co-occurrence control._
  * GL_ ripples_hfos _control.m

_Mentioned in Text: Shuffling Plusmaze co-occurrence control._
  * GL_plusmaze_control.m

_**Figure 3 (A,B,C,D)**: Spectral power during events._
  * GL_spectral_power.m  (Older version 2020)
  * GL_spectral_power_2021.m (Updated version 2022)

_**Figure 3 (E,F,G)**: Granger causality during events._
  * GL_granger.m
  * GL_granger_Nayanika.m (Computes spectral Granger Causality analysis after combining events of all rats).

_**Figure 4A**:_ 
  * GL_delta_counts.m

_**Figure 4B**: Spindles counts_
  * GL_spindles_counts.m * (Older version 2020)
  * GL_spindles_counts_nayanika.m * (Updated version 2022)

_**Figure 4C**:_ 
  * GL_delta_spindles.* (Computes the co-occurrence of deltas an spindles from PPC and PFC as done by [Kim et al.,Cell, 2019](https://www.cell.com/cell/pdf/S0092-8674(19)30959-6.pdf))

_**Figure 4 (D,E,F)**: Spindle co-occurrence. Before & After counts._	
  * GL_spindles.m * (Older version 2020)
  * GL_spindles_Nayanika.m * (Updated version 2022)

_Mentioned in Text: Spindle co-occurrence shuffling control_
  * GL_spindles_control.m *
 
_**Figure 4I**:_ 
  * GL_swr_disruption.m
 
   
---------
*In the current version (2022): Spindles were computed using an adaptation of the FindSpindles.m function from [FMA toolbox](https://github.com/michael-zugaro/FMAToolbox/tree/master/Analyses) (Zugaro lab).

*In the older version (2020): Spindles were previously detected using the YASA algorithm. The steps to do this were:
1. Run GL_spindle_matlab2python.m for every session per condition to export NREM epochs to python.
2. Run GL_yasa_spindles.py for every session per condition to save detections in a .mat file.

--------------------------------

## Outdated fork.
[Plusmaze 2022 update](https://github.com/nayanikab20/GL_fast_slow_hfos) (Nayanika) 
