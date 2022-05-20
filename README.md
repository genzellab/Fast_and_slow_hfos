# Fast and slow cortical high frequency oscillations for cortico-cortical and corticohippocampal network consolidation during NonREM sleep. 

Adrian Aleman-Zapata 1, Richard GM Morris 2, Lisa Genzel *1,2

*corresponding author: lgenzel@donders.ru.nl  :mailbox: 

1, Donders Institute for Brain Cognition and Behavior, Radboud University, Postbus 9010, 6500GL Nijmegen/Netherlands.

2, Centre for Cognitive and Neural Systems, Edinburgh Neuroscience, University of Edinburgh, 1 George Square, Edinburgh EH8 9JZ, UK.

Reference:  [doi.org/10.1101/765149](https://doi.org/10.1101/765149) 

-----------------------------




:warning: Requirements: Makes use of functions from the [YASA](https://github.com/raphaelvallat/yasa) , [Fieldtrip](https://github.com/fieldtrip/fieldtrip) and [ADRITOOLS](https://github.com/Aleman-Z/ADRITOOLS) repositories. 

These last two need to be added to the path.

We also include in the subfunctions folder some functions taken from [FMA toolbox](https://github.com/michael-zugaro/FMAToolbox/tree/master/Analyses)

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

_Figure 1 and 2: Count of coocurring and single events._
  * GL_hfos_counts.m

_Figure 2C: Shuffling co-occurrence control._
  * GL_ ripples_hfos _control.m

_Figure 2D: Shuffling Plusmaze co-occurrence control._
  * GL_plusmaze_control.m

_Figure 3: Spectral power during events._
  * GL_spectral_power.m 
  * GL_spectral_power_2021.m

_Figure 4A: Spindles counts_
  * GL_spindles_counts.m *

_Figure 4 (B,C,D,G): Spindle co-occurrence. Before & After counts._	
  * GL_spindles.m *

_Figure 4 (F,G): Spindle co-occurrence shuffling control_
  * GL_spindles_control.m *

_Figure 6: Granger causality during events._
  * GL_granger.m
  
_New (2021):_ 
  * GL_swr_disruption.m

_Newer (2022):_ 
  * GL_delta_counts.m
  * GL_spindles_counts_nayanika.m *
  * GL_spindles_Nayanika.m *
  * GL_delta_spindles.* (Computes the co-occurrence of deltas an spindles from PPC and PFC as done by [Kim et al.,Cell, 2019](https://www.cell.com/cell/pdf/S0092-8674(19)30959-6.pdf))
  * GL_granger_Nayanika.m (Computes spectral Granger Causality analysis after combining events of all rats).
  
---------
*In the current version (2022): Spindles were computed using an adaptation of the FindSpindles.m function from [FMA toolbox](https://github.com/michael-zugaro/FMAToolbox/tree/master/Analyses) (Zugaro lab).

*In the older version (2020): Spindles were previously detected using the YASA algorithm. The steps to do this were:
1. Run GL_spindle_matlab2python.m for every session per condition to export NREM epochs to python.
2. Run GL_yasa_spindles.py for every session per condition to save detections in a .mat file.

--------------------------------
