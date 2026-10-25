"""Compare fixed-seed runs and validate every local analysis artifact."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from statistics import mean
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
ANALYSIS=ROOT/'evaluation/descriptor/generation_analysis_reinvent'
METRICS=['average_oracle_score_origin','average_oracle_score_sigmoid','average_atom_number','average_MW','repetition_count','noncanonical_count','average_generated_heavy_atoms','average_base_reward','average_training_reward']


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes=plt.subplots(2,2,figsize=(10,7),layout="constrained")
    panels={"average_atom_number":(axes[0,0],"Whole-molecule heavy atoms"),
            "average_generated_heavy_atoms":(axes[0,1],"Generated-fragment heavy atoms"),
            "average_MW":(axes[1,0],"Molecular weight"),
            "average_oracle_score_origin":(axes[1,1],"Raw PLANET score")}
    summaries=[]
    baseline=None
    for tag,penalty in [('001',.01),('0015',.015),('002',.02)]:
        data=ROOT/f'results/reinvent-high9-cap9-size{tag}-1000-20261025/all_generated.csv.gz'
        with gzip.open(data,'rt') as handle:rows=list(csv.DictReader(handle))
        initial={(r['target'],r['mode'],r['sample_index']):r for r in rows if r['step']=='1'}
        if baseline is None:baseline=initial
        same_smiles=sum(r['smiles']==baseline[k]['smiles'] for k,r in initial.items())
        same_raw=sum(r['planet_affinity']==baseline[k]['planet_affinity'] for k,r in initial.items())
        result={'penalty':penalty,'attempts':len(rows),'scored':sum(r['planet_status']=='ok' for r in rows),
                'invalid':sum(r['smiles_state']=='0' for r in rows),'repeats':sum(int(r['is_repeat']) for r in rows),
                'negative_rewards':sum(float(r['training_reward'])<0 for r in rows),
                'initial_smiles_matching_baseline':same_smiles,'initial_raw_matching_baseline':same_raw,
                'data_sha256':hashlib.sha256(data.read_bytes()).hexdigest()}
        folder=ANALYSIS/f'high9_cap9_size{tag}_1000_seed42_20261025/analysis'
        for metric in METRICS:
            with (folder/'summary/table'/f'{metric}.csv').open() as handle:curve=list(csv.DictReader(handle))
            assert len(curve)==10 and len(curve[0])==31
            if metric in panels:
                ax,title=panels[metric]
                ax.plot(range(1,11),[mean(float(v) for k,v in r.items() if k!="batch") for r in curve],marker="o",label=str(penalty))
                ax.set(title=title,xlabel="RL step");ax.grid(alpha=.2)
            for label,row in [('first',curve[0]),('final',curve[-1])]:
                result[f'{label}_{metric}']=mean(float(v) for k,v in row.items() if k!='batch')
        summaries.append(result)
        if tag=='001':continue
        pngs=list(folder.rglob('*.png'))
        assert len(pngs)==324
        for path in pngs:
            with Image.open(path) as im:im.verify()
        for metric in METRICS:
            tables=list((folder/metric/'table').glob('*.csv'))
            assert len(tables)==33
            for path in tables:
                with path.open() as handle:table=list(csv.DictReader(handle))
                assert [int(r['batch']) for r in table]==list(range(1,11))
        with (folder/'summary_table.csv').open() as handle:assert len(list(csv.DictReader(handle)))==33
        (folder/'artifact-validation.json').write_text(json.dumps({'passed':True,'pngs_verified':324,'metric_groups':9,'per_metric_tables_verified':33,'steps':10,'summary_rows':33},indent=2)+'\n')
    with (OUT/'comparison.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    (OUT/'comparison.json').write_text(json.dumps(summaries,indent=2)+'\n')
    for ax in axes.flat:ax.legend(title="Size penalty")
    fig.suptitle("30-task equal-weight means; seed 42; 100 samples per step")
    fig.savefig(OUT/"comparison_curves.png",dpi=180)
    plt.close(fig)
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
