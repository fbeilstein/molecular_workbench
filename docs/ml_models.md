# Machine Learning Transition State Models (RitS & TS-Diff)

Molecular Workbench supports optional neural-network generative transition state predictors such as **RitS** (Reaction-informed Transition States via Flow-Matching / Diffusion) from the Isayev Lab / MIT.

---

## Why Isn't RitS Pre-Installed?

1. **Heavyweight Dependencies**: Machine learning flow-matching models require PyTorch (`torch`), PyTorch Geometric (`torch-geometric`), and specific CUDA/C++ extensions (~2.5 GB download).
2. **Model Weights Checkpoints**: Trained neural network weights checkpoints (`rits_weights.pt` or `transition1x.pt`) are 1.2 to 2.0 GB in size and are hosted on external model hubs (Zenodo / Hugging Face).
3. **No Silent Fallbacks**: The workbench will **never** silently fake a neural network calculation if the package or weights are missing. If selected without prerequisites, it raises a clear error pointing to this setup guide.

---

## How to Install RitS

### 1. Install PyTorch & PyG into your environment:
```bash
# In your active .venv:
pip install torch --index-url https://download.pytorch.org/whl/cpu  # or cu121 for GPU
pip install torch-geometric
```

### 2. Install the RitS package:
Run `setup.sh --with-ml` or install directly from `tools/rits_repo`:
```bash
pip install --no-deps -e tools/rits_repo
```

### 3. Download Model Weights Checkpoint:
The official RitS model checkpoints are hosted by the authors on Google Drive:
**Official Release Folder**: [Google Drive Release Folder](https://drive.google.com/drive/folders/1xGX0uZgaGHU_e-a8hDdflV4wWApczNVU?usp=sharing)

Files available:
- `rits.ckpt`: Checkpoint trained on the full RGD1-ChEMBL-extended dataset (~1.7M reactions).
- `ts1x_rits.ckpt`: Checkpoint trained on the Transition1x benchmark.

#### Option A: Download via browser
1. Open the [Google Drive link](https://drive.google.com/drive/folders/1xGX0uZgaGHU_e-a8hDdflV4wWApczNVU?usp=sharing) and download `rits.ckpt`.
2. Move it to:
   ```bash
   tools/ml_models/rits/rits.ckpt
   ```

#### Option B: Download via `gdown`
```bash
pip install gdown
mkdir -p tools/ml_models/rits
gdown "https://drive.google.com/drive/folders/1xGX0uZgaGHU_e-a8hDdflV4wWApczNVU" -O tools/ml_models/rits/ --folder
```

### 4. Automatic Discovery:
The engine automatically detects checkpoints at any of these paths:
- `tools/ml_models/rits/rits.ckpt`
- `tools/ml_models/rits/rits_checkpoint.pt`
- `tools/rits_repo/data/rits.ckpt`
- Or via environment variable: `export RITS_CHECKPOINT="/path/to/rits.ckpt"`
