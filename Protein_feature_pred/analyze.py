import os
import re
import subprocess
from pathlib import Path

import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_minced_and_prodigal(input_dir, output_dir):
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for fna_file in input_dir.glob("*.fna"):
        basename = fna_file.stem

        output_txt = output_dir / f"{basename}.txt"
        output_gff = output_dir / f"{basename}.gff"
        output_faa = output_dir / f"{basename}.faa"
        output_gene_gff = output_dir / f"{basename}_gene.gff"

        cmd_minced = ["./minced", str(fna_file), str(output_txt), str(output_gff)]
        print("Running minced:", " ".join(cmd_minced))
        subprocess.run(cmd_minced)

        cmd_prodigal = [
            "prodigal", "-i", str(fna_file), "-c", "-m", "-f", "gff",
            "-a", str(output_faa), "-o", str(output_gene_gff), "-p", "meta"
        ]
        print("Running prodigal:", " ".join(cmd_prodigal))
        subprocess.run(cmd_prodigal)

def parse_minced_gff(gff_file):
    data = []
    basename = Path(gff_file).stem

    with open(gff_file) as handle:
        for record in GFF.parse(handle):
            for feature in record.features:
                if feature.type == "repeat_region":
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    rpt_seq = feature.qualifiers.get("rpt_unit_seq", [""])[0]
                    n_repeat = feature.qualifiers.get("score", [""])[0]
                    data.append({
                        "MAG": basename,
                        "Contig": record.id,
                        "CRISPR_Start": start,
                        "CRISPR_End": end,
                        "Direct repeats": rpt_seq,
                        "N_repeats": n_repeat,
                        "Spacers": n_repeat
                    })
    return pd.DataFrame(data)

def parse_minced_txt(txt_file):
    basename = Path(txt_file).stem
    # 读取文件内容
    with open(txt_file, 'r') as file:
        raw_text = file.read()

    # 分割不同序列
    sequences = re.split(r"Sequence '", raw_text)
    parsed_data = {}

    for seq in sequences:
        if not seq.strip():
            continue

        lines = seq.strip().splitlines()
        seq_id_match = re.match(r"([^']+)' \((\d+) bp\)", lines[0])
        if not seq_id_match:
            continue

        seq_id, length = seq_id_match.groups()
        current_crispr = None
        data = []

        for line in lines[1:]:
            if m := re.match(r"CRISPR (\d+)\s+Range: (\d+) - (\d+)", line):
                current_crispr = {
                    "crispr_id": int(m[1]),
                    "range": (int(m[2]), int(m[3])),
                    "repeats": []
                }
            elif m := re.match(r"(\d+)\s+([ATCG]+)\s+([ATCG]+)\s+\[\s*(\d+),\s*(\d+)\s*\]", line):
                position, repeat, spacer, rlen, slen = m.groups()
                current_crispr["repeats"].append({
                    "position": int(position),
                    "repeat": repeat,
                    "spacer": spacer,
                    "repeat_len": int(rlen),
                    "spacer_len": int(slen)
                })
            elif m := re.match(r"(\d+)\s+([ATCG]+)", line):
                position, repeat = m.groups()
                current_crispr["repeats"].append({
                    "position": int(position),
                    "repeat": repeat,
                    "spacer": "",
                    "repeat_len": len(repeat),
                    "spacer_len": 0
                })
            elif line.startswith("Repeats:") and current_crispr:
                parsed_data[f"{seq_id}"] = current_crispr

    data=[]
    # 打印结构化数据
    for key, value in parsed_data.items():
        Contig = key
        CRISPR_Start = value['range'][0]
        CRISPR_End = value['range'][1]
        Direct_Repeats = value['repeats'][0]["repeat"]
        N_Repeats = len(value['repeats'])
        N_Spacers = N_Repeats
        
        Forward_chain_start = value['repeats'][0]["position"]+value['repeats'][0]["repeat_len"]+value['repeats'][0]["spacer_len"]
        Reverse_chain_start = value['repeats'][-1]["position"]

        if value['repeats'][-1]["spacer_len"] == 0:
            N_Spacers=N_Spacers-1
            Reverse_chain_start = Reverse_chain_start - value['repeats'][-2]["spacer_len"]
        data.append({"MAG":basename,"Contig":Contig,"CRISPR_Start":CRISPR_Start,"CRISPR_End":CRISPR_End,"Direct_Repeats":Direct_Repeats,"N_Repeats":N_Repeats,"N_Spacers":N_Spacers,"Forward_chain_start":Forward_chain_start,"Reverse_chain_start":Reverse_chain_start})
        '''
        print(f"CRISPR标识符: {key}")
        print(f"  编号: {value['crispr_id']}")
        print(f"  区间: {value['range']}")
        print("  重复单元:")
        for repeat_info in value['repeats']:
            print(f"    位置: {repeat_info['position']}")
            print(f"    REPEAT: {repeat_info['repeat']}")
            print(f"    SPACER: {repeat_info['spacer']}")
            print(f"    REPEAT长度: {repeat_info['repeat_len']}, SPACER长度: {repeat_info['spacer_len']}")
        print()
        '''
    # df = pd.DataFrame(data)
    # print(df)
    return pd.DataFrame(data)

def parse_protein_fasta(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        start, end,pn_chain = record.description.split("#")[1:4]
        data.append({
            "id": record.id,
            "start": start,
            "end": end,
            "pn_chain":pn_chain,
            "seq": str(record.seq)
        })
    return pd.DataFrame(data)

def generate_candidate_fasta(minced_txt_file, protein_faa_file, output_fasta_file):
    df_crispr = parse_minced_txt(minced_txt_file)
    df_proteins = parse_protein_fasta(protein_faa_file)

    contigs = df_crispr['Contig'].unique()
    filtered = df_proteins[df_proteins['id'].apply(lambda x: any(x.startswith(c) for c in contigs))]

    records = [
        SeqRecord(Seq(row["seq"]), id=row["id"], description=f"{row['start']}-{row['end']}")
        for _, row in filtered.iterrows()
    ]

    SeqIO.write(records, output_fasta_file, "fasta")

def predit_castype_with_model():
    '''
    accelerate launch --mixed_precision=fp16 \
    --use_deepspeed --config_file config/deepspeed_config.yaml \
    esm2_classification.py -a predict \
    --eval_dataset_dir datasets/candidate \
    -m epoch_0/ \
    -l 1560 --label_file config/labels.txt \
    --output_file candidate.csv
    '''

def get_subtype(keys_str,stat):

    present_keys = [key.strip() for key in keys_str.split('_') if key.strip()]
    if not present_keys:
        return "输入不能为空"
    
    # 检查必须存在的Cas12
    if 'Cas12' not in present_keys:
        return "other"
    
    idx = present_keys.index("Cas12")
    if stat>0:
        present_keys = present_keys[:idx + 1]
    else:
        present_keys = present_keys[idx:]

    
    # 定义所有可能的组件（Cas1, Cas2, Cas4, Cas12）
    all_components = {'Cas1', 'Cas2', 'Cas4', 'Cas12'}
    # 计算各组件的存在状态（存在为1，不存在为0）
    cas_status = {
        'Cas1': 1 if 'Cas1' in present_keys else 0,
        'Cas2': 1 if 'Cas2' in present_keys else 0,
        'Cas4': 1 if 'Cas4' in present_keys else 0,
        'Cas12': 1  # 已验证Cas12存在
    }
    
    # 提取状态值
    cas1, cas2, cas4 = cas_status['Cas1'], cas_status['Cas2'], cas_status['Cas4']
    
    # 根据规则分类
    if cas1 and cas2 and cas4:
        return "1"
    elif not cas1 and cas2 and cas4:
        return "2"
    elif cas1 and not cas2 and cas4:
        return "3"
    elif cas1 and cas2 and not cas4:
        return "4"
    elif not cas1 and not cas2 and cas4:
        return "5"
    elif cas1 and not cas2 and not cas4:
        return "6"
    elif not cas1 and cas2 and not cas4:
        return "7"
    elif not cas1 and not cas2 and not cas4:
        return "8"
    else:
        return "未定义的组合"  # 理论上不会出现，因为状态由三个二进制位完全枚举


def merge_prediction_results(prediction_file, minced_txt_file, protein_faa_file, output_file):
    df_pred = pd.read_csv(prediction_file)
    df_pred = df_pred[df_pred["predicted_label"] != "nocas"]

    df_crispr = parse_minced_txt(minced_txt_file)
    df_proteins = parse_protein_fasta(protein_faa_file)

    df_merged = pd.merge(df_pred, df_proteins, left_on="name", right_on="id", how="inner")

    results = []
    for _, row in df_crispr.iterrows():
        matching = df_merged[df_merged["name"].str.startswith(row["Contig"])]
        print(matching)
        cas12_end = matching[matching["predicted_label"]=="Cas12"]
        
        # loci length计算为Cas12起始到第一个spacer结束

        print(cas12_end)
        if not cas12_end.empty:
            operon_start = int(matching.iloc[0]["start"])
            operon_end = int(cas12_end.iloc[-1]["end"])
            genes = "_".join(matching["predicted_label"])
            proteins = ",".join(matching["name"])

            loci_length = -1
            forward_chain_start = row["Forward_chain_start"]
            reverse_chain_start = row["Reverse_chain_start"]
            pn_chain = matching.iloc[0]["pn_chain"]
            if pn_chain == -1:
                loci_length = abs(operon_end-reverse_chain_start)+1
            else:
                loci_length = abs(operon_end-forward_chain_start)+1

            results.append({
                **row,
                "Operon_Start": operon_start,
                "Operon_End": operon_end,
                "CRISPR-loci type":get_subtype(genes,operon_start-row['CRISPR_End']),
                "CRISPR-loci length":loci_length,
                "Genes": genes,
                "Proteins": proteins,
            })

    pd.DataFrame(results).to_csv(output_file, index=False)

# Example Usage (uncomment and adapt as needed):
# run_minced_and_prodigal("input", "output")
# generate_candidate_fasta("output/GMBC10.002_143.txt", "output/GMBC10.002_143.faa", "output/candidate.faa")
# predit_castype_with_model()
merge_prediction_results("output/candidate.csv", "output/GMBC10.002_143.txt", "output/GMBC10.002_143.faa", "output/final_summary.csv")
# parse_minced("output/info.txt")
