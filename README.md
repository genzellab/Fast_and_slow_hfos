# Fast and slow cortical high frequency oscillations for cortico-cortical and corticohippocampal network consolidation during NonREM sleep. 

Adrian Aleman-Zapata 1, Richard GM Morris 2, Lisa Genzel *1,2

*corresponding author: lgenzel@donders.ru.nl  :mailbox: 

1, Donders Institute for Brain Cognition and Behavior, Radboud University, Postbus 9010, 6500GL Nijmegen/Netherlands.

2, Centre for Cognitive and Neural Systems, Edinburgh Neuroscience, University of Edinburgh, 1 George Square, Edinburgh EH8 9JZ, UK.

Reference:  [doi.org/10.1101/765149](https://doi.org/10.1101/765149) 

-----------------------------




:warning: Requirements: Makes use of functions from the [YASA](https://github.com/raphaelvallat/yasa) , [Fieldtrip](https://github.com/fieldtrip/fieldtrip) and [ADRITOOLS](https://github.com/Aleman-Z/ADRITOOLS) repositories. 
These last two need to be added to the path.

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

_Figure 4A: Spindles counts_
  * GL_spindles_counts.m *

_Figure 4 (B,C,D,G): Spindle co-occurrence. Before & After counts._	
  * GL_spindles.m *

_Figure 4 (F,G): Spindle co-occurrence shuffling control_
  * GL_spindles_control.m *

_Figure 6: Granger causality during events._
  * GL_granger.m


---------
*Spindles were previously detected using the YASA algorithm. The steps to do this were:
1. Run GL_spindle_matlab2python.m for every session per condition to export NREM epochs to python.
2. Run GL_yasa_spindles.py for every session per condition to save detections in a .mat file.

--------------------------------
