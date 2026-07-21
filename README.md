# Community Fast-SL 🧬

A computationally efficient Python package for identifying synthetic lethality in microbial communities using MICOM and the Fast-SL pruning algorithm.

## 🔬 Overview
This tool performs exhaustive *in silico* Double Lethal (DL) knockouts across two different species in a shared environment. It drastically reduces execution time by mathematically pruning "Zero-Flux $\times$ Zero-Flux" combinations using Parsimonious FBA (pFBA), avoiding hundreds of thousands of redundant Gurobi solver iterations.

## 🚀 Installation
Clone this repository and install the package in editable mode:

```bash
git clone [https://github.com/yash-yadav-official/Community_FastSL.git](https://github.com/yash-yadav-official/Community_FastSL.git)
cd Community_FastSL
pip install -e .
```
## 📖 Interactive Tutorial
We have provided a complete, step-by-step interactive tutorial for running an exhaustive screen between E. coli and Salmonella.
You can access the jupyter notebook tutorial by running ```
jupiter notebook```in the terminal in the community_fastSL directory.

you can also just click the tutorial notebook in github to look into it. furthermore you can open the folder in vscode and use from there as well. 
Use the tutorial as a template to work on the fastSL_workbench notebook.
