# EpPredictor
Enhancer-Promoter interactions(EPI) predictions using machine learning approach 

Welcome to the EpPredictor, Here is the file format:

- core---save the codes for extract epigenomic featrues

  - muti_feature_extractor.py---extract feature from chip-seqs and save results
  - utils.py---tools functions

- data--the folder store chip-seq data

  - GM12878
  - HUVEC
  - ...
  - protine_info.xlsx---detail of chip-seqs we used

- EP-interaction---store sequence data, extract sequence features and model, visualization..etc

  - data---save extracted chip-seq features, sequence features and EP interaction pairs.

    - EP-lower---store extracted sequence features

    - GM12878

      - midfile--store interaction pairs
      - muti-midfile---store extracted features from various region(window, promoter2k)
      - train.csv---completed features for training model

    - HUVEC

      ...

    - sequence---store source sequence data

      - hg19.id_chr10.fa
      - hg19.id_chr1.fa
      - ...

  - scripts---codes for extract sequence features

    - ep_lower.py---sequence feature extractor
    - merge_all.py---merge sequence feature and chip-seq features
    - run_eplower.sh---bash for run ep_lower.py
    - tools.py---some help functions

  - Visualization---some jupyter scripts for train model and draw pictures

## Train the model based on extracted features

Prepare: The train files should be stored in **path/EpPredictor/EP-interaction/cell line/muti-midfile/train.csv**make sure the train files has been stored correctly, you can obtain these files from:

[requires data](https://drive.google.com/open?id=1Wih5l07BnQ47r6kQUtCywFKaaEiRi-K3), download the data and replace folder **path/EpPredictor/EP-interaction/data**

then into folder **path/EpPredictor/EP-interaction/Viaualization/model/**, and open **VariousModel2.ipynb** via jupyter notebook, this file contain various models use various region feature to predict EP interactions.

## Generate features and train model step by step

- step1: prepare souce data, including interaction pairs, chip-seqs and sequences

  1. all the chip-seq feature can be download from ENCODE Project or NCBI, the detail chip-seq feature we used can be find in **protein_info.xlsx**, it's not necessary use complete chip-seqs, you can use as much as you can get. Next, put these files into *path/EpPredictor/data/cell line/xxx.bed*
  2. download interaction pairs from [interactions](https://drive.google.com/open?id=1Wih5l07BnQ47r6kQUtCywFKaaEiRi-K3) or [targetfinder](https://github.com/shwhalen/targetfinder.git), and put these file into *path/EpPredictor/EP-interactions/data/cell line/midfile/xxx.csv*
  3. download sequence data(fastaq format), put it into *path/EpPredictor/EP-interaction/data/sequence/hg19.id_chrxx.fa*

- step2: extract chip-seq features:

  into *path/EpPredictor/core/*, and run **muti_feature_extractor.py**:

  Usage: ```python muti_feature_extractor.py --chipseq_path '../data/IMR90/ --interaction_path '../EP-interaction/data/IMR90/midfile/IMR90_pairs.csv' --feature_save_path   '../EP-interaction/data/IMR90/muti-midfile/'```

  | parameter name    |      | default                                                | meaning                 |
  | ----------------- | ---- | ------------------------------------------------------ | ----------------------- |
  | chipseq_path      |      | '../data/IMR90/                                        | chip-seq data folder    |
  | interaction_path  |      | '../EP-interaction/data/IMR90/midfile/IMR90_pairs.csv' | interaction file        |
  | feature_save_path |      | '../EP-interaction/data/IMR90/muti-midfile/'           | where to store features |

  

  make sure all the parameters including file path was right. This step could be taken times, when finished, you can obtain the chip-seq features in feature paths.

- step3: merge chip seq features:

  into folder *path/EpPredictor/EP-interaction/scripts/*, and run **merge_all.py** 
 : ```python run merge_all.py ```

- step4: extract sequence features

  into folder *path/EpPredictor/EP-interaction/scripts*, and run command: 

  **sh run_eplower.sh**

- step5: all features has been extracted, you can train model now like above

  

  









