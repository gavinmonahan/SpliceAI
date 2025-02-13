## Instructions to setup tensorflow-gpu on Setonix

Following (this guide)[https://pawsey.atlassian.net/wiki/spaces/US/pages/51924596/Running+TensorFlow+on+Setonix#RunningTensorFlowonSetonix-InstallingadditionalPythonpackages]

1. Disable conda
`nano ~/.bashrc`
2. Swap modules
```
module unload gcc/12.2.0
module swap pawseyenv/2024.05 pawseyenv/2023.08
module load gcc/12.2.0
module load tensorflow/rocm5.6-tf2.12
```
3. bash into Singulrity container
```
mkdir -p $MYSOFTWARE/manual/software/pythonEnvironments/tensorflowContainer-environments
cd $MYSOFTWARE/manual/software/pythonEnvironments/tensorflowContainer-environments
bash
cd $MYSOFTWARE/manual/software/pythonEnvironments/tensorflowContainer-environments
```
4. Create venv
```
python3 -m venv --system-site-packages spliceai2
python3 -m pip install /software/projects/pawsey0933/gmonahan/T2T/SpliceAI/
```
5. Check installation
`/home/gmonahan/.local/bin/spliceai`
