# 🧠 Fast and Slow Cortical High-Frequency Oscillations During NonREM Sleep

This repository contains code and analyses for the project:

**Fast and slow cortical high frequency oscillations for cortico-cortical and cortico-hippocampal network consolidation during NonREM sleep.**  
📖 [Preprint on bioRxiv](https://doi.org/10.1101/765149)

This project later became the published paper:

**Sleep deprivation and hippocampal ripple disruption after one-session learning eliminate memory expression the next day.**  
📖 [PNAS paper](https://doi.org/10.1073/pnas.2123424119)



**Authors:**  
Adrián Alemán-Zapata¹, Richard G.M. Morris², Lisa Genzel*¹²  
\*Corresponding author: [lgenzel@donders.ru.nl](mailto:lgenzel@donders.ru.nl)

**Affiliations:**  
1. Donders Institute for Brain, Cognition and Behaviour, Radboud University, Nijmegen, Netherlands  
2. Centre for Cognitive and Neural Systems, University of Edinburgh, UK

---

## ⚙️ Requirements

This repository relies on:

### 📦 Toolboxes and Dependencies

- **MATLAB toolboxes:**
  - Image Processing Toolbox
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Mapping Toolbox
  - Deep Learning Toolbox
  - Symbolic Math Toolbox
  - Bioinformatics Toolbox
  - Computer Vision Toolbox
  - Fixed-Point Designer
  - MATLAB Coder
  - Simulink
  - Parallel Computing Toolbox
  - MATLAB Parallel Server
  - Polyspace Bug Finder

- **External libraries:**
  - [FieldTrip](https://github.com/fieldtrip/fieldtrip) — tested with 2019 version to avoid plotting issues
  - [ADRITOOLS](https://github.com/Aleman-Z/ADRITOOLS)
  - [FMA Toolbox](https://github.com/michael-zugaro/FMAToolbox) — specific analysis functions included in `subfunctions/`
  - [YASA](https://github.com/raphaelvallat/yasa) (used in the 2020 version only)

> ⚠️ FieldTrip and ADRITOOLS must be added to the MATLAB path manually.

---

## 📂 Main Scripts by Figure

### 🧮 Co-occurrence and Event Counting

- **Figure 2B, 2D** – Count co-occurring/single events, slow/fast counts and rates:  
  `GL_hfos_counts.m`

- **Text Analysis – Ripple-HFO shuffling controls:**  
  `GL_ripples_hfos_control.m`  
  `GL_plusmaze_control.m`

---

### 📈 Spectral Power & Granger Causality (Figure 3)

- **Figure 3A–D** – Spectral power:  
  `GL_spectral_power.m` (2020 version)  
  `GL_spectral_power_2021.m` (2022 version)

- **Figure 3E–G** – Granger causality:  
  `GL_granger.m`  
  `GL_granger_Nayanika.m` (all rats combined)

---

### 🌊 Delta and Spindle Analyses (Figure 4)

- **Figure 4A** – Delta counts:  
  `GL_delta_counts.m`

- **Figure 4B** – Spindle counts:  
  `GL_spindles_counts.m` (2020)  
  `GL_spindles_counts_nayanika.m` (2022)

- **Figure 4C** – Delta–Spindle co-occurrence (inspired by [Kim et al., *Cell*, 2019](https://www.cell.com/cell/pdf/S0092-8674(19)30959-6.pdf)):  
  `GL_delta_spindles.m`

- **Figure 4D–F** – Spindle co-occurrence (before/after counts):  
  `GL_spindles.m` (2020)  
  `GL_spindles_Nayanika.m` (2022)

- **Text Analysis – Spindle co-occurrence control:**  
  `GL_spindles_control.m`

- **Figure 4I** – SWR disruption:  
  `GL_swr_disruption.m`

---

## 🔄 Spindle Detection Versions

- **2022 (Current):**  
  Spindles detected using a modified version of `FindSpindles.m` from the Zugaro Lab's [FMA Toolbox](https://github.com/michael-zugaro/FMAToolbox).

- **2020 (Old):**  
  Spindles detected using YASA (Python):
  1. Export NREM epochs using `GL_spindle_matlab2python.m`
  2. Run spindle detection with `GL_yasa_spindles.py`  
     Results are saved as `.mat` files

---

## 🧪 Related Work (Outdated)

You can find an outdated fork focused on Plusmaze analysis here:  
🔗 [Plusmaze 2022 fork by Nayanika](https://github.com/nayanikab20/GL_fast_slow_hfos)

---

## 📬 Contact

For questions or collaboration inquiries, feel free to reach out:

- 👨‍🔬 Adrian Aleman-Zapata: [GitHub](https://github.com/Aleman-Z) | [Email](mailto:adrian.alemanzapata@donders.ru.nl)
- 📧 Lisa Genzel: [lgenzel@donders.ru.nl](mailto:lgenzel@donders.ru.nl)

---

