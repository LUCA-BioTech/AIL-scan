
# Overview

Cas classification and Cas12a trans-cleavage activity prediction

## Cas classification:
This project focuses on the classification of CRISPR-Cas system proteins. We obtained Cas protein sequences from the NCBI database and fine-tuned the ESM model to build an efficient prediction tool for accurately classifying different types of Cas proteins. This model can be utilized for in-depth studies of CRISPR system mechanisms and serves as a reference for the development of gene-editing tools.

## Cas12a trans-cleavage activity prediction:
Cas12a trans-cleavage activity prediction and visualization
Small sample learning methods and results for Cas12a trans-cleavage activity prediction and visualization.
Users need to set a path to save the visualization results.
Three dimensionality reduction methods can be used for comparison, including 'umap', 'pca' and 'tsne'.
Train and test data are provided, but only contain embeddings output by our fine-tuned ESM model for confidentiality and patent reasons.

### Important entrances:
* file esm2_classification.py : Cas classification code
* forder trans_active : Cas12a trans-cleavage activity prediction
* forder compare : Cas classification compare , support esm nakh prot_t5
* forder Cas12 structure : 3D structure of Cas12 candidate proteins
* forder protein-embeddings : Attention visualize and TNE code

## Requirements
```
conda env create -f AIL-scan.yml

The num_processes in config/deepspeed_config.yml should match the number of available GPUs.
```

## Cas classification General usage

### Fetching data
Run `collect_dataset.py` to fetch the CRISPR-associated protein dataset. 

### train / eval / predict / export

```
$ python esm2_classification.py

optional arguments:
    * parser.add_argument('--train_dataset_dir', type=str, help='train datasets directory that includes FASTA files named according to their labels')
    * parser.add_argument('--eval_dataset_dir', type=str, required=True, help='eval datasets directory that includes FASTA files named according to their labels')
    * parser.add_argument('--label_file', type=str, help='label file for train/eval/predict')
    * parser.add_argument('-a', '--action', type=str, required=True, help='action', choices=["train", "eval", "predict", "export"])
    * parser.add_argument('--attention_strategy', type=str, help='calculation strategy for attention export', default=None, choices=["full-huge-storage-to-layer", "full-huge-storage-to-head", "2d", "3d-max", "3d-average", "3d-high_confidence_average"])
    * parser.add_argument('--attention_layers', nargs='+', type=str, help='layers for attention export, use "all" for all layers', default=[-3, -2, -1])
    * parser.add_argument('--attention_threshold', type=float, help='threshold for attention calculation')
    * parser.add_argument('--predict_head_mask', type=parse_head_mask, help='Mask to nullify selected heads of the self-attention modules. Mask values selected in [0, 1]: 1 indicates the head is not masked, 0 indicates the head is masked. examples: --predict_head_mask  "[0, 1, ... all heads ..., 1, 0]"  or --predict_head_mask "[[0, 1, ... all heads ..., 1, 0], ... all layers ..., [1, 0, ... all heads ..., 0, 1]]"')
    * parser.add_argument('-m', '--model_name', type=str, required=True, help='model name')
    * parser.add_argument('-l', '--max_seq_len', type=int, required=True, help='max sequence length')
    * parser.add_argument('-o', '--output_dir', type=str, help='output directory')
    * parser.add_argument('--num_train_epochs', type=int, help='number of epochs', default=1)
    * parser.add_argument('--learning_rate', type=float, help='learning rate', default=1e-5)
    * parser.add_argument('--min_learning_rate', type=float, help='min learning rate', default=0)
    * parser.add_argument('--max_train_steps', type=int, help='max_train_steps', default=None)
    * parser.add_argument('--micro_train_batch_size', type=int, help='train batch size per device', default=2)
    * parser.add_argument('--eval_batch_size', type=int, help='eval batch size per device', default=10)
    * parser.add_argument('--random_seed', type=int, help='random seed', default=42)
    * parser.add_argument('--gradient_checkpointing', action='store_true', help='gradient checkpointing')
    * parser.add_argument('--gradient_accumulation_steps', type=int, help='gradient accumulation steps', default=1)
    * parser.add_argument('--weight_decay', type=float, help='weight_decay', default=0.01)
    * parser.add_argument('--save_steps', type=int, help='save_steps', default=-1)
    * parser.add_argument('--num_warmup_steps', type=int, help='warmup steps', default=1000)
    * parser.add_argument('--logging_steps', type=int, help='logging_steps', default=100)
    * parser.add_argument('--logging_dir', type=str, help='logging_dir', default="./logs")
    * parser.add_argument('--output_file', type=str, help='predict/export output file')
    * parser.add_argument('--checkpoint', type=str, help="checkpoint for resume")
    * parser.add_argument('--loss_fn', type=str, help="Loss function", choices=["focal"])
    * parser.add_argument("--project", type=str, help="W&B project name", default="esm2_classification")
    * parser.add_argument('--local_rank', type=int, help='local rank', default=-1)
```

For example:

### train a model
```
accelerate launch --mixed_precision=fp16 \
  --use_deepspeed --config_file config/deepspeed_config.yaml \
  esm2_classification.py -a train --train_dataset_dir datasets/train \
  --eval_dataset_dir datasets/validate \
  -m /home/luca/.cache/huggingface/hub/models--facebook--esm2_t33_650M_UR50D \
  -l 1560 \
  --micro_train_batch_size 4 \
  -o ./output \
  --label_file config/labels.txt \
  --num_train_epochs 5 --num_warmup_steps 500 --weight_decay 0.0001 \
  --learning_rate 0.001
```

### predict
```
accelerate launch --mixed_precision=fp16 \
 --use_deepspeed --config_file config/deepspeed_config.yaml \
 esm2_classification.py -a predict \
 --eval_dataset_dir datasets/test \
 -m output/models--facebook--esm2_t33_650M_UR50D/2024-12-10_15-08-08/epoch_0 \
 -l 1560 --label_file config/labels.txt \
 --output_file 2024-12-10_15-08-08-epoch_0.csv
```

### export attention
```
accelerate launch --mixed_precision=fp16 \
 --use_deepspeed --config_file config/deepspeed_config.yaml \
 esm2_classification.py -a export \
 --eval_dataset_dir datasets/export \
 -m output/models--facebook--esm2_t33_650M_UR50D/2024-12-02_00-47-40/epoch_0/ \
 -l 1560 --label_file config/labels.txt \
 --output_dir export \
 --eval_batch_size 1 \
 --micro_train_batch_size 1 \
 --attention_strategy full-huge-storage-to-layer \
 --attention_threshold 0.05
```

##  Cas12a trans-cleavage activity prediction General usage

Small sample learning methods and results for Cas12a trans-cleavage activity prediction and visualization.

### Code for training: pipeline.py
Users need to set a path to save the visualization results. (The default save_path is None)
Three dimensionality reduction methods can be used for comparison, including 'umap', 'pca' and 'tsne'. 
Train and test data are provided, but only contain embeddings output by our fine-tuned ESM for confidentiality and patent reasons.
Models will be saved in folder 'AutoglounModels'; pca pickle files will be saved in folder 'AutoglounRd'; results will be saved in folder 'AutoglounLogs'.
The running time is around 158 seconds for training process with 2 4090 GPUs. (Running pipeline with training and testing data and save all models for 2-15 dimensions.)
Directly run pipeline.py under the required environment.

### Code for testing: test.py 
We have prepared a model with one of the best performances along with the code and test.py will automatically read models and pca.pkl file in folder autogloun and pca respectively.
We provide 2 test sets for testing. One is newtest.csv which contains the two protein mentioned in our paper. (GMBC10.001_783__k119_63415_106 and GMBC10.008_184__k119_67936_25) The other is newtest1.csv which contains 14 protein sequences. The embeddings are output from the fine-tuned ESM.
The running time is around 61 seconds for predicting trans-cleavage activity for newtest1.csv (14 protein sequences) with saved models and 2 4090 GPUs.
Directly run test.py under the required environment. 0 indicates no trans-cleavage activity and 1 indicates with trans-cleavage activity. Change the "new_test_path" for different testing data and tasks.

# Third-Party Dependencies

This project incorporates the following open-source tools and libraries. Their licenses and copyright notices are detailed in the `LICENSES` directory:

| Tool/Library Name | License Type | Copyright Notice | Code Path |
|----------------------|-------------------|--------------------------|---------------------|
| **minced** | [GPL v3](https://github.com/ctSkennerton/minced?tab=GPL-3.0-1-ov-file) | Copyright (c) 2011-2019  Florent Angly, Connor Skennerton | [GitHub](https://github.com/ctSkennerton/minced) |
| **prodigal** | [GPL v3](https://github.com/hyattpd/Prodigal?tab=GPL-3.0-1-ov-file) | Copyright (C) 2010-2023 Doug Hyatt et al. | [GitHub](https://github.com/hyattpd/Prodigal) |
| **esm2** | [MIT](https://github.com/facebookresearch/esm?tab=MIT-1-ov-file) | Copyright (c) 2021 Facebook AI Research (FAIR) | [GitHub](https://github.com/facebookresearch/esm) |
| **deepspeed** | [Apache 2.0](https://github.com/deepspeedai/DeepSpeed?tab=Apache-2.0-1-ov-file)  | Copyright (c) 2020 Microsoft Corporation | [GitHub](https://github.com/deepspeedai/DeepSpeed) |
| **accelerate** | [Apache 2.0](https://github.com/huggingface/accelerate?tab=Apache-2.0-1-ov-file) | Copyright (c) 2022 Hugging Face Inc. | [GitHub](https://github.com/huggingface/accelerate) |


## License Compatibility Statement

This project is licensed under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), which is compatible with the licenses of all dependencies listed above:

- **prodigal** (GPL v3) is used as an independent tool without modifying its source code. According to GPL v3 terms, this does not require changing the license of this project.
- Other dependencies (MIT/Apache 2.0) permit integration into Apache 2.0 projects.


## Usage Recommendations

If you modify or distribute this project, ensure compliance with the licenses of all dependencies and retain the above copyright notices and license texts.


# License
This source code is licensed under the Apache-2.0 license found in the LICENSE file in the root directory of this source tree.





