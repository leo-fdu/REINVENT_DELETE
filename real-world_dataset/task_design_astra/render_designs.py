#!/usr/bin/env python3
"""Reproduce the design atlas; never write model inputs or modify raw ligands.

Reference structures are atom-mapped to MOL2 atom IDs. RCSB CCD templates
(downloaded 2026-09-22) supply bond orders, hydrogens and stereochemistry for
15 topology-matched ligands. DRD3 retains the local graph and uses its 3D stereo.
Full-ligand reference strings below are provenance, not model input records.
"""
from pathlib import Path
import hashlib
import html
import io
import sys

import torch
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT / "REINVENT4"))
from reinvent.chemistry import conversions
from reinvent.chemistry.library_design import attachment_points, bond_maker
from reinvent.models.transformer.core.vocabulary import SMILESTokenizer

CODES = {'aa2ar': 'ZMA',
 'abl1': 'JIN',
 'aces': 'HUX',
 'adrb1': 'P32',
 'akt1': '0XZ',
 'cdk2': 'FAP',
 'cp3a4': 'RIT',
 'dpp4': 'KIQ',
 'drd3': 'ETQ',
 'dyr': 'D2B',
 'egfr': 'HYZ',
 'gria2': 'ZK1',
 'jak2': 'NVB',
 'kif11': 'K30',
 'parp1': 'A92',
 'ppara': '735'}

REFERENCES = {'aa2ar': '[cH:1]1[cH:2][c:3]([OH:4])[cH:5][cH:6][c:7]1[CH2:8][CH2:9][NH:10][c:11]1[n:12][c:18]2[n:16]([c:14]([NH2:15])[n:13]1)[n:17][c:20](-[c:21]1[cH:22][cH:23][cH:24][o:25]1)[n:19]2',
 'abl1': '[CH3:1][c:11]1[cH:8][c:12]([NH:24][c:21]2[n:22][cH:10][c:16]3[cH:9][c:17](-[c:18]4[c:14]([Cl:28])[cH:4][cH:3][cH:5][c:15]4[Cl:29])[c:20](=[O:26])[n:25]([CH3:2])[c:19]3[n:23]2)[cH:6][cH:7][c:13]1[F:27]',
 'aces': '[CH3:1][CH2:9][C:18]1=[CH:7][C@@H:16]2[CH2:5][c:13]3[c:14]([c:15]([NH2:19])[c:11]4[cH:3][cH:2][c:10]([Cl:21])[cH:4][c:12]4[n:20]3)[C@@H:17]([CH2:6]2)[CH2:8]1',
 'adrb1': '[OH:1][C@@H:2]([CH2:4][NH:5][C:6]([CH3:7])([CH3:8])[CH3:21])[CH2:9][O:10][c:11]1[cH:12][cH:13][cH:14][c:15]2[c:16]1[CH2:19][C:18]([C:20]#[N:3])=[N:17]2',
 'akt1': '[Cl:1][c:2]1[cH:3][cH:4][c:7]([C@H:8]([CH2:9][CH2:10][OH:11])[NH:12][C:13](=[O:14])[C:15]2([NH2:18])[CH2:16][CH2:17][N:21]([c:22]3[n:23][cH:24][n:25][c:30]4[c:26]3[cH:27][cH:28][nH:29]4)[CH2:20][CH2:19]2)[cH:6][cH:5]1',
 'cdk2': '[CH3:1][N:26]([CH3:2])[CH2:12][C@@H:21]([CH2:13][O:28][c:15]1[cH:8][cH:6][c:14]([NH:24][c:18]2[cH:10][c:19]([NH:25][c:20]3[c:16]([F:29])[cH:4][cH:3][cH:5][c:17]3[F:30])[n:23][cH:11][n:22]2)[cH:7][cH:9]1)[OH:27]',
 'cp3a4': '[CH3:1][CH:29]([CH3:2])[c:28]1[n:39][c:27]([CH2:21][N:43]([CH3:5])[C:36]([NH:41][C@@H:33]([CH:30]([CH3:3])[CH3:4])[C:35]([NH:42][C@@H:34]([CH2:20][c:25]2[cH:14][cH:10][cH:7][cH:11][cH:15]2)[CH2:23][C@@H:32]([C@H:31]([CH2:19][c:24]2[cH:12][cH:8][cH:6][cH:9][cH:13]2)[NH:40][C:37](=[O:47])[O:48][CH2:22][c:26]2[cH:16][n:38][cH:18][s:50]2)[OH:44])=[O:45])=[O:46])[cH:17][s:49]1',
 'dpp4': '[cH:1]1[c:9]([C@H:15]2[CH2:5][CH:8]=[C:17]([C:18]([N:24]3[CH2:4][CH2:3][n:23]4[c:13]([n:21][n:22][c:14]4[C:19]([F:29])([F:30])[F:31])[CH2:6]3)=[O:25])[CH2:7][C@@H:16]2[NH2:20])[c:11]([F:27])[cH:2][c:12]([F:28])[c:10]1[F:26]',
 'drd3': '[CH3:1][c:10]1[cH:8][c:5]([CH2:22][CH3:23])[c:2]([OH:4])[c:13]([C:14]([NH:3][CH2:15][C@H:16]2[N:6]([CH2:20][CH3:21])[CH2:19][CH2:18][CH2:17]2)=[O:9])[c:11]1[O:7][CH3:12]',
 'dyr': '[CH3:1][CH:17]([CH3:2])/[C:18](=[CH:9]/[c:12]1[cH:8][o:23][c:15]2[c:13]1[c:14]([NH2:19])[n:21][c:16]([NH2:20])[n:22]2)[c:10]1[cH:6][cH:4][cH:5][cH:7][c:11]1[O:24][CH3:3]',
 'egfr': '[cH:1]1[cH:2][c:17]([CH2:16][n:31]2[c:20]3[cH:5][cH:4][c:18]([NH:30][c:24]4[c:22]([CH2:15][NH:29][N:32]5[CH2:13][CH2:11][CH2:10][CH2:12][CH2:14]5)[c:23]([NH2:25])[n:27][cH:9][n:28]4)[cH:7][c:21]3[cH:8][n:26]2)[cH:6][c:19]([F:33])[cH:3]1',
 'gria2': '[cH:1]1[c:8]([C:14]([F:24])([F:25])[F:26])[c:11]([N:17]2[CH2:3][CH2:5][O:23][CH2:6][CH2:4]2)[cH:2][c:10]2[c:9]1[nH:15][c:12](=[O:20])[c:13](=[O:21])[n:16]2[CH2:7][P:27](=[O:18])([OH:19])[OH:22]',
 'jak2': '[CH3:1][NH:27][S:33]([c:16]1[cH:10][cH:8][c:15](-[c:17]2[cH:6][cH:5][cH:7][c:18]3[c:23]2[n:26][c:22](-[c:19]2[cH:12][c:20]([O:30][CH3:2])[c:24]([O:32][CH3:4])[c:21]([O:31][CH3:3])[cH:13]2)[cH:14][n:25]3)[cH:9][cH:11]1)(=[O:28])=[O:29]',
 'kif11': '[CH3:1][N:26]1[CH2:12][CH2:11][C@H:21]([N:28]([CH3:2])[C:24]([N:27]2[CH2:15][C:22]([c:20]3[cH:10][c:18]([F:31])[cH:8][cH:9][c:19]3[F:32])=[CH:13][C@@:25]2([CH2:16][OH:29])[c:17]2[cH:6][cH:4][cH:3][cH:5][cH:7]2)=[O:30])[C@H:23]([F:33])[CH2:14]1',
 'parp1': '[cH:1]1[cH:2][c:13]([C:19]([NH2:20])=[O:24])[c:16]2[c:14]([cH:3]1)[n:22][c:17](-[c:12]1[cH:5][cH:4][c:11]([C@@H:18]3[CH2:9][CH2:7][CH2:8][CH2:10][NH:23]3)[cH:6][c:15]1[F:25])[nH:21]2',
 'ppara': '[CH3:1][c:13]1[c:18]([C:20]([NH:25][CH2:12][c:14]2[cH:4][cH:10][c:17]([O:29][C:22]([CH3:2])([CH3:3])[C:21](=[O:26])[OH:28])[cH:11][cH:5]2)=[O:27])[s:33][c:19](-[c:15]2[cH:6][cH:8][c:16]([C:23]([F:30])([F:31])[F:32])[cH:9][cH:7]2)[n:24]1'}

RAW_SHA256 = {'aa2ar': 'fb4f772045e80d15e919fb0cc748f4e13429e73d9b6d284ab203f51233042cc7',
 'abl1': '351685615a31721c552492a3af427e4e5ac55f50061a67a41917457403efc0e6',
 'aces': '1703c50c57a2c24592b2e6fb571be4b7aa3bf1604e8d1ebd591b231cbfb2631e',
 'adrb1': 'cf21c52fd56a44a837fc25989a185f07d3880e3e158e6a9e08a80c622ad0d0c7',
 'akt1': '4ad33e8a4308897c232f5b08213a1917a17b29a307c5e7fa6827338a9e6098ae',
 'cdk2': '721c460ebf2db564aae3b85953a5548971b476fa8e0c0863adfe1c2c67c57c6a',
 'cp3a4': 'af38b4da66f376d216054a89849690ad6d90d6150d6e05ee99b7c65f11099bdb',
 'dpp4': 'c66e797c433aa245a2ac2af30527e0f12178f38c08dcbe1aecf247d564b0a843',
 'drd3': 'd6168ae6539f9de16c471c1dc2696757d6dd5269c65bfd122650aa93f80e2594',
 'dyr': '16136044517330b22775c8c859d1a39ecf431d97b5bb2f5173ab99d9eb68a142',
 'egfr': 'c19cce10b32333f6a4a03b7e072865adc1c40ab5c83f8ad8bdc27882c3863bf8',
 'gria2': '1ef46e9db70370659bad4a1a32c36ca6f14594cc719653791534f54837a2cea1',
 'jak2': 'cf672194d8e746a28a01968b8b30124b89423340d3b9c61e7ddfdad9dc3a80a6',
 'kif11': '70a9d35ee57ce2752f7e707d63102b5eefc9fb916ce59a871b5d2ddcd316b8fa',
 'parp1': '495f20545fafcd06d76c0edf88aec52d0ec85337ab6806ea54ee5e9d57143e8a',
 'ppara': '35db7f1d001b3d5723e2817f576394ac8436044a6890708689d2ebaedc29fd7a'}

# MOL2 atom IDs: LibINVENT cut, retained-side seed; LinkINVENT two cuts.
DESIGNS = {
 'aa2ar': ((10,9),11,[(10,9),(8,7)], '保留含氨基的稠合杂环和呋喃，重建苯乙基侧链。', '保留杂环端与苯酚端，替换二碳连接链。'),
 'abl1': ((12,24),21,[(12,24),(24,21)], '保留含二氯苯基的杂环核心，装饰苯胺端。', '保留两侧环系，替换单原子 NH 桥；属于短 linker 任务。'),
 'aces': ((18,9),18,[(18,9),(9,1)], '保留完整多环骨架，在原乙基位置生长侧链。', '边界案例：保留多环核心与甲基，仅替换 CH2。甲基仅 1 个重原子，不建议纳入主 linking 基准。'),
 'adrb1': ((5,6),11,[(10,9),(4,5)], '保留含氰基环系和氨基醇链，重新装饰胺端。', '保留含氰基环系的醚端与叔丁胺端，替换羟基丙链。'),
 'akt1': ((13,12),15,[(13,12),(7,8)], '保留吡咯并嘧啶和氨基哌啶酰基，重建酰胺侧链。', '保留酰基核心与氯苯端，替换带羟乙基支链的含氮连接区。'),
 'cdk2': ((28,13),15,[(14,24),(18,24)], '保留二芳胺嘧啶和芳氧基，生长含氨基醇的侧链。', '保留两侧环系，替换单原子 NH 桥；属于短 linker 任务。'),
 'cp3a4': ((37,48),31,[(31,32),(34,23)], '保留主链及另一噻唑端，装饰末端氨基甲酸酯的氧侧。', '保留两侧较大片段，替换 CH(OH)-CH2 连接区。'),
 'dpp4': ((18,24),17,[(17,18),(18,24)], '保留氨基环己烯、氟苯和酰基，重建含氮杂环端。', '保留氨基环己烯芳基端与含氮杂环端，替换羰基短桥，避免在手性中心切断。'),
 'drd3': ((3,15),13,[(3,14),(16,15)], '按本地配体保留取代苯甲酰胺端，重建吡咯烷侧链。', '按本地配体保留芳酰基和 N-乙基吡咯烷，替换 NH-CH2。'),
 'dyr': ((24,11),12,[(9,12),(18,10)], '保留二氨基杂环、支化烯链和苯基，重建甲氧基取代位，避免切断双键邻位。', '保留二氨基杂环与邻甲氧基苯基，替换支化烯链；当前拼接会丢失参考 E/Z 信息。'),
 'egfr': ((31,16),24,[(17,16),(31,16)], '保留含氮多环与氨基嘧啶区域，装饰 N-苄基端。', '保留含氮核心和氟苯端，替换单个 CH2 桥；属于短 linker 任务。'),
 'gria2': ((17,11),16,[(16,7),(27,7)], '保留二酮稠环和膦酸侧链，在芳碳位重建吗啉端。', '保留稠环与膦酸基，替换单个 CH2 桥；膦酸端较小，需单独分层。'),
 'jak2': ((19,22),17,[(15,17),(19,22)], '保留含磺酰胺的芳基与稠合杂环，装饰三甲氧基苯端。', '保留两端取代苯基，替换中间稠合杂环；这是含环 linker / 核心替换任务。'),
 'kif11': ((20,22),25,[(27,24),(21,28)], '保留含氟哌啶、脲基与主杂环，装饰二氟苯端。', '保留两侧环系，替换 C(=O)-N(CH3) 连接区。'),
 'parp1': ((11,18),17,[(12,17),(11,18)], '保留苯并咪唑甲酰胺和氟苯，重建哌啶端。', '保留苯并咪唑甲酰胺与哌啶端，替换中间氟苯连接区。'),
 'ppara': ((19,15),18,[(18,20),(12,14)], '保留噻唑、酰胺与羧酸端，装饰三氟甲基苯端。', '保留噻唑芳基端和含羧酸的苯氧基端，替换 C(=O)-NH-CH2。'),
}
PDB = dict(zip(CODES, '3EML 2HZI 1E66 2VT4 4GV1 1H00 3NXU 2I78 3PBL 3NXO 2RGP 3KGC 3LPB 3CJO 3L3M 2P54'.split()))
COLORS = [(0.65,0.85,0.98),(0.69,0.89,0.74),(1.0,0.80,0.54)]
FONT = '/System/Library/Fonts/Supplemental/Arial Unicode.ttf'

def load_ligand(target):
    path = ROOT / 'real-world_dataset' / target / 'crystal.mol2'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == RAW_SHA256[target], path
    mol = Chem.MolFromSmiles(REFERENCES[target])
    order = sorted(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetAtomMapNum())
    mol = Chem.RenumberAtoms(mol, order)
    for atom in mol.GetAtoms():
        assert atom.GetAtomMapNum() == atom.GetIdx() + 1
        atom.SetIntProp('source_id', atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    # Verify the archived reference still represents this local heavy-atom graph.
    lines = path.read_text().splitlines()
    atom_start = lines.index('@<TRIPOS>ATOM') + 1
    raw_atoms = []
    for line in lines[atom_start:]:
        if line.startswith('@<TRIPOS>'):
            break
        fields = line.split()
        raw_atoms.append((int(fields[0]), fields[5].split('.')[0]))
    assert raw_atoms == [(a.GetIdx()+1,a.GetSymbol()) for a in mol.GetAtoms()]
    bond_start = lines.index('@<TRIPOS>BOND') + 1
    raw_edges = set()
    for line in lines[bond_start:]:
        if line.startswith('@<TRIPOS>'):
            break
        if line.strip():
            fields = line.split()
            raw_edges.add(tuple(sorted((int(fields[1]),int(fields[2])))))
    assert raw_edges == {tuple(sorted((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1))) for b in mol.GetBonds()}
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    rdDepictor.Compute2DCoords(mol)
    return mol

def split(mol, cuts):
    bonds = [mol.GetBondBetweenAtoms(a-1,b-1) for a,b in cuts]
    assert all(b is not None and not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE for b in bonds)
    fragmented = Chem.FragmentOnBonds(mol, [b.GetIdx() for b in bonds], dummyLabels=[(i+1,i+1) for i in range(len(cuts))])
    fragments = list(Chem.GetMolFrags(fragmented, asMols=True))
    for f in fragments:
        for a in f.GetAtoms():
            if a.GetAtomicNum() == 0:
                a.SetIntProp('cut_id', a.GetIsotope())
                a.SetIsotope(0)
        Chem.SanitizeMol(f)
        assert not any(a.GetNumRadicalElectrons() for a in f.GetAtoms())
    return fragments

def ids(mol):
    return {a.GetIntProp('source_id') for a in mol.GetAtoms() if a.GetAtomicNum()}

def dummies(mol):
    return [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]

def smi(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def clean(mol):
    mol = Chem.Mol(mol)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    return smi(mol)

def check_vocab(text, model, label):
    tokens = SMILESTokenizer().tokenize(text)
    unknown = sorted(set(tokens) - set(model['vocabulary'].tokens()))
    assert not unknown, (label, unknown)
    assert len(tokens) <= model['max_sequence_length'], (label,len(tokens))
    return len(tokens)

def validate(mol, core, removed, ends, linker, models, target):
    assert len(dummies(core)) == len(dummies(removed)) == 1
    assert len(ends) == 2 and all(len(dummies(e)) == 1 for e in ends)
    assert len(dummies(linker)) == 2
    for f in [core,removed,*ends,linker]:
        assert len(Chem.GetMolFrags(f)) == 1
        assert all(a.GetDegree() == 1 and a.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE for a in dummies(f))
        assert Chem.MolFromSmiles(smi(f)) is not None
    # Use the repository's actual assembler; reference fragments only in memory.
    numbered_core = attachment_points.add_attachment_point_numbers(smi(core), canonicalize=False)
    lib_joined = bond_maker.join_scaffolds_and_decorations(numbered_core,smi(removed),keep_labels_on_atoms=True)
    assert lib_joined is not None and clean(lib_joined) == clean(mol), 'LibINVENT reconstruction'
    # Preserve the relationship between the linker's SMILES traversal and warhead order.
    labeled_linker = Chem.Mol(linker)
    for a in dummies(labeled_linker): a.SetAtomMapNum(a.GetIntProp('cut_id'))
    parsed = Chem.MolFromSmiles(smi(labeled_linker))
    cut_order = [a.GetAtomMapNum() for a in dummies(parsed)]
    ordered = [next(e for e in ends if dummies(e)[0].GetIntProp('cut_id') == cut) for cut in cut_order]
    for a in parsed.GetAtoms(): a.SetAtomMapNum(0)
    linker_text = Chem.MolToSmiles(parsed, canonical=False, isomericSmiles=True)
    joined = bond_maker.join_scaffolds_and_decorations(
        attachment_points.add_attachment_point_numbers(linker_text,canonicalize=False),
        '|'.join(smi(e) for e in ordered), keep_labels_on_atoms=True)
    assert joined is not None, 'LinkINVENT assembly failed'
    if target == 'dyr':
        # Existing assembler drops E/Z at this cut. Require exact connectivity,
        # and explicitly report the missing stereochemistry in the figure/docs.
        no_maps = Chem.Mol(joined)
        for a in no_maps.GetAtoms(): a.SetAtomMapNum(0)
        assert Chem.MolToSmiles(no_maps,isomericSmiles=False) == Chem.MolToSmiles(mol,isomericSmiles=False)
        assert clean(joined) != clean(mol), 'Update the documented E/Z limitation if assembler changes'
    else:
        assert clean(joined) == clean(mol), 'LinkINVENT stereochemical reconstruction'
    # Standardization follows Transformer sampling; no model input file is emitted.
    lib_input = conversions.convert_to_standardized_smiles(smi(core))
    link_input = '|'.join(conversions.convert_to_standardized_smiles(smi(e)) for e in ends)
    for original, standardized in [(core,lib_input),*zip(ends,link_input.split('|'))]:
        assert clean(Chem.MolFromSmiles(standardized)) == clean(original), 'Standardization changed input'
    token_counts = [check_vocab(lib_input,models['libinvent'],'libinvent'),check_vocab(link_input,models['linkinvent'],'linkinvent')]
    check_vocab(smi(removed), models['libinvent'], 'reference decoration')
    check_vocab(linker_text, models['linkinvent'], 'reference linker')
    return token_counts

def font(size):
    return ImageFont.truetype(FONT,size)

def text(draw, xy, message, size=25, color='#243447', width=70):
    # Wrap Chinese by measured pixel width; width denotes approximate Latin chars.
    line=''; y=xy[1]
    for ch in message:
        if ch=='\n' or draw.textlength(line+ch,font=font(size)) > width*size*0.52:
            draw.text((xy[0],y),line,font=font(size),fill=color); y+=size+9; line=''
        if ch!='\n': line+=ch
    draw.text((xy[0],y),line,font=font(size),fill=color)
    return y+size+9

def draw_molecule(mol, size, color_by_id=None, cuts=(), labels=False):
    mol=Chem.Mol(mol)
    if not color_by_id:
        rdDepictor.Compute2DCoords(mol)
    drawer=rdMolDraw2D.MolDraw2DCairo(*size)
    opt=drawer.drawOptions(); opt.padding=0.10; opt.minFontSize=15; opt.maxFontSize=28
    atom_colors={}; bond_colors={}
    if color_by_id:
        atom_colors={a.GetIdx():color_by_id[a.GetIntProp('source_id')] for a in mol.GetAtoms()}
        for b in mol.GetBonds():
            a,z=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
            if atom_colors[a]==atom_colors[z]:bond_colors[b.GetIdx()]=atom_colors[a]
    for a,b in cuts:
        bond_colors[mol.GetBondBetweenAtoms(a-1,b-1).GetIdx()] = (0.94,0.36,0.37)
        for n in [a,b]:mol.GetAtomWithIdx(n-1).SetProp('atomNote',str(n))
    if labels:
        for a in dummies(mol): opt.atomLabels[a.GetIdx()]='*'
    drawer.DrawMolecule(mol,highlightAtoms=list(atom_colors),highlightBonds=list(bond_colors),highlightAtomColors=atom_colors,highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    return Image.open(io.BytesIO(drawer.GetDrawingText())).convert('RGB')

def render(target,mol,core,removed,ends,linker,counts):
    cut,seed,linkcuts,lib_reason,link_reason=DESIGNS[target]
    canvas=Image.new('RGB',(2100,1640),'#f1f5f9'); draw=ImageDraw.Draw(canvas)
    text(draw,(48,22),f'{target.upper()}  |  配体任务设计',44,width=90)
    source=f'本地 crystal.mol2 · {PDB[target]} / {CODES[target]} · 重原子 {mol.GetNumHeavyAtoms()}'
    if target=='drd3':source='本地 crystal.mol2 · 保留本地结构（与 3PBL / ETQ 的 Cl→CH3 差异）'
    text(draw,(50,83),source,24,width=135)
    text(draw,(50,121),'蓝 / 绿：固定保留片段    橙：待生成区域    红：切断键    *：模型连接点；邻近数字为 MOL2 原子 ID',24,width=155)
    for row,kind,frags,cutlist,reason in [(0,'LibINVENT',[core,removed],[cut],lib_reason),(1,'LinkINVENT',[*ends,linker],linkcuts,link_reason)]:
        top=177+row*690
        draw.rounded_rectangle((28,top,2072,top+668),radius=20,fill='white')
        text(draw,(52,top+17),kind,34,width=90)
        text(draw,(305,top+22),reason,24,width=130)
        colors={}
        for i,f in enumerate(frags):
            col=COLORS[2] if i==len(frags)-1 else COLORS[i]
            colors.update({j:col for j in ids(f)})
        text(draw,(57,top+103),'参考配体：切分示意',24,width=60)
        canvas.paste(draw_molecule(mol,(880,420),colors,cutlist),(45,top+141))
        draw.line((948,top+101,948,top+594),fill='#d9e2ec',width=3)
        if row==0:
            panels=[(core,'固定骨架 · 1 个 *',970,580),(removed,'参考装饰片段 · 待生成',1575,465)]
        else:
            panels=[(ends[0],'固定片段 A · 1 个 *',970,350),(ends[1],'固定片段 B · 1 个 *',1340,350),(linker,'参考 linker · 待生成',1710,330)]
        for f,title,x,w in panels:
            text(draw,(x+8,top+105),title,22,width=32)
            canvas.paste(draw_molecule(f,(w,400),labels=True),(x,top+147))
            text(draw,(x+10,top+550),f'{f.GetNumHeavyAtoms()} 个重原子',21,width=30)
        if row==1:
            text(draw,(1305,top+331),'|',38,color='#51677e',width=3)
        cutlabel='；'.join(f'{a}—{b}' for a,b in cutlist)
        status='输入连接点 / 词表 / 参考拼回检查通过'
        if target=='dyr' and row==1:status='输入检查通过；拼回仅连接图一致，E/Z 丢失'
        if target=='aces' and row==1:status='形式检查通过；甲基端过小，仅作边界示例'
        text(draw,(57,top+610),f'切断键（MOL2 原子 ID）：{cutlabel}    |    {status}    |    输入 token 数 {counts[row]}/180',22,width=170)
    text(draw,(49,1571),'仅展示设计，不是模型输入文件。橙色结构是原配体中的参考答案，不作为生成条件；二维切分不代表已验证结合贡献。',22,width=175)
    canvas.save(OUT/f'{target}.png',optimize=True)

def main():
    models={k:torch.load(ROOT/f'REINVENT4/priors/{k}_transformer_pubchem.prior',map_location='cpu',weights_only=False) for k in ['libinvent','linkinvent']}
    rows=[]
    for target in CODES:
        mol=load_ligand(target)
        cut,seed,linkcuts,lib_reason,link_reason=DESIGNS[target]
        parts=split(mol,[cut]);assert len(parts)==2
        core=next(f for f in parts if seed in ids(f));removed=next(f for f in parts if f is not core)
        parts=split(mol,linkcuts);assert len(parts)==3
        ends=[f for f in parts if len(dummies(f))==1]
        linker=next(f for f in parts if len(dummies(f))==2)
        assert set.union(*(ids(f) for f in parts))==ids(mol)
        assert sum(len(ids(f)) for f in parts)==len(ids(mol))
        counts=validate(mol,core,removed,ends,linker,models,target)
        render(target,mol,core,removed,ends,linker,counts)
        rows.append((target,core.GetNumHeavyAtoms(),removed.GetNumHeavyAtoms(),*[e.GetNumHeavyAtoms() for e in ends],linker.GetNumHeavyAtoms(),counts))
        limitation = '; LinkINVENT E/Z lost on reconstruction' if target == 'dyr' else ''
        print(f'{target}: validated and drawn; heavy atoms {rows[-1][1:6]}; tokens {counts}{limitation}')
    make_docs(rows,models)


def make_docs(rows,models):
    links=' '.join(f'<a href="#{t}">{t.upper()}</a>' for t in CODES)
    cards='\n'.join(f'<section id="{t}"><h2>{t.upper()}</h2><p>{html.escape(DESIGNS[t][3])}<br>{html.escape(DESIGNS[t][4])}</p><a href="{t}.png"><img loading="lazy" src="{t}.png" alt="{t} LibINVENT 和 LinkINVENT 二维任务设计"></a></section>' for t in CODES)
    (OUT/'index.html').write_text('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>16 靶点 · LibINVENT / LinkINVENT 任务设计</title><style>body{margin:0;background:#edf2f7;color:#223449;font:17px/1.6 system-ui,sans-serif}header,main{max-width:1450px;margin:auto;padding:24px}header{background:white}nav{display:flex;flex-wrap:wrap;gap:16px}a{color:#155b9a}section{background:white;padding:22px;margin:24px 0;border-radius:12px}img{width:100%;height:auto}h1{font-size:30px}.note{border-left:4px solid #de962e;padding:12px;background:#fff5e4}</style><header><h1>16 靶点 · LibINVENT / LinkINVENT 二维任务设计</h1><p>每个靶点各 1 个 scaffold decoration / growing 方案和 1 个 fragment linking 方案。蓝绿为固定条件，橙色为待生成区域，红色为切断键。点击图片可查看原尺寸。</p><p class="note">ACES linking 是单重原子端片段的边界案例，不建议纳入主基准；短 linker 与含环 linker 应分层比较。DRD3 保留本地配体，与 PDB 标准配体存在元素差异。原 MOL2 的键级与氢表示问题已在绘图参考结构中核对，原文件未改动。</p><p><a href="README.md">设计依据、数据差异和验证说明</a></p><nav>'''+links+'</nav></header><main>'+cards+'</main></html>',encoding='utf-8')
    table='\n'.join(f'| [{t}]({t}.png) | {a} + {b} | {c} + {d} + {e} | {n[0]} / {n[1]} |' for t,a,b,c,d,e,n in rows)
    sources='\n'.join(f'- `{t}`：[{PDB[t]}](https://www.rcsb.org/structure/{PDB[t]}) / [{CODES[t]}](https://www.rcsb.org/ligand/{CODES[t]})；原始 SHA-256 `{RAW_SHA256[t]}`。' for t in CODES)
    revisions='\n'.join(f'- `{k}`：metadata model_id `{v["metadata"]["model_id"]}`；max_sequence_length `{v["max_sequence_length"]}`。' for k,v in models.items())
    readme=f'''# 16 个靶点的二维任务设计

打开 [index.html](index.html) 浏览全部图；也可直接查看各靶点 PNG。每张图上半为 LibINVENT，下半为 LinkINVENT，共 32 个设计。这里只交付图、说明和复现绘图脚本，没有生成 `.smi`、SDF/MOL2 片段或运行配置。

## 输入规则与图例

- **LibINVENT**：一个连续骨架，在固定的出口放一个 `*`；本设计统一使用单出口，便于与 DELETE growing 比较。模型生成另一侧装饰片段，不将原装饰片段作为输入。
- **LinkINVENT**：两个独立且各带一个 `*` 的片段。实际 REINVENT4 记录用 `片段A|片段B` 的形式；不是点号分隔的 SMILES。图中竖线表达该分隔符。生成 linker 有两个出口，不能把一个连续骨架上的两个出口当作这两个输入片段。
- 蓝色/绿色为保留原子，橙色为原配体中待替换区域（参考答案）；红色表示切断键。所有切断均为非环单键，不破坏芳香环或切断双键。`*` 是 dummy attachment atom，不是真实碳原子、不计重原子数。
- 图中连接点不加化学同位素。切断键使用原 MOL2 的 **1-based atom ID**，图上在对应端点显示这些 ID。LibINVENT 内部重新编号、LinkINVENT 的片段顺序与 linker 出口配对由实际 sampler 处理，不能把图上的原子 ID 直接当作模型的连接点编号。

## 任务价值与比较范围

这些是按结构提出的初始任务，不是经蛋白接触分析确认的最佳药效团切分。没有运行生成、打分、对接或亲和力验证。每个靶点的理由见图与图册。

- **ACES**：保留完整多环骨架和末端甲基，替换乙基中的一个 CH2。模型格式与拼接验证通过，但甲基只有 1 个重原子，不能视为有充分结合贡献的独立片段。保留为边界案例，建议主 linking 汇总排除它并另报覆盖率；若强制要求两个较大的片段，应另选已知配体，而不是在此多环骨架上硬切环。
- **ABL1、CDK2、DPP4、EGFR、GRIA2**：参考 linker 为单个 NH、CH2 或羰基（GRIA2 的另一端是较小膦酸基）。格式有效，单独报告短 linker 层；不能与长链重建任务直接混为一种难度。
- **JAK2、PARP1**：参考 linker 含环，属于含环连接区/核心替换任务。LinkINVENT 支持其连接拓扑，但生成难度和固定药效团范围与柔性链替换不同，应单列结果。
- **DYR linking**：输入满足连接点和词表要求，参考答案可恢复全部键连接；但当前 REINVENT4 拼接实现丢失该切口对应的烯键 E/Z 信息。它不是立体化学完全还原的任务，需另加几何异构体约束或单列结果。LibINVENT 方案已选用远离烯键的切口，保留参考 E/Z。
- 其余任务同样记录实际保留和删除重原子数。跨模型比较必须让 DELETE 与相应 REINVENT4 任务使用完全相同的原子保留集合；DELETE 的片段坐标将来应从原晶体坐标提取，不能使用二维绘图坐标。
- 原配体恢复仅作完整性检查；正式基准应另外报告新颖性，避免把恢复已知配体当作创新。模型输入的二维结构并不自行施加蛋白结合位置或三维出口方向约束。

## 化学结构核对

原始 MOL2 没有显式氢，RDKit 直接读取有 7 个无法 sanitize，其他若干会出现异常自由基、价态或电荷。不能使用 `sanitize=False` 得到的图就宣称可作为模型输入。

本图册对 15 个配体采用 RCSB CCD 字典的完整配体参考结构，在**元素和重原子邻接图完全一致**的前提下映射回原 MOL2 原子 ID；用字典的键级、氢、形式电荷和立体化学绘图。CCD 理想坐标没有用于任务片段空间约束。参考字符串及原文件哈希嵌入脚本，复现不需要联网。原始数据未修改。

这不是对原 MOL2 进行无条件认可：例如 CP3A4、DYR、KIF11 等的键级与字典有差异；ADRB1 字典使用其记录中的亚胺/非芳香五元环互变异构表示；羧酸、膦酸与胺使用参考结构的电荷态。这里采用可解析的参考状态，并不等同于某一 pH 的主导质子化态。以后真正构建输入前，应确认两模型统一采用同一个化学标准化方案。

**DRD3 特例**：本地文件原子 1 是连接芳环的 C，3PBL 的 ETQ 字典相应位置为 Cl，元素图不能完整匹配。没有用 ETQ 替换本地配体。本图保留本地全部 23 个重原子的元素、键级和连接，补足正常价态的隐式氢，胺采用中性表示，手性取自本地三维坐标。图的任务设计对本地结构成立，但不可把它直接称作标准 ETQ 共晶配体。

## 验证

绘图脚本对全部 32 个方案执行：

1. 原 MOL2 SHA-256 与本次审阅版本一致；所有参考结构可 sanitize，无异常自由基。
2. 非环单键切分；LibINVENT 得到 2 个连续片段，固定端 1 个出口；LinkINVENT 得到 3 个连续片段，两个固定端各 1 个出口，待生成区 2 个出口。
3. 原子无遗漏或重复，dummy 均为度 1、单键连接。
4. 使用本仓库 `bond_maker.join_scaffolds_and_decorations` 和连接点编号逻辑，将参考答案拼回，比较去除 atom-map 后的 canonical isomeric SMILES；31 个方案完全一致，DYR linking 仅连接图一致、E/Z 丢失，已单独标注。单原子 linker 也经过实际拼接检查。
5. 使用 Transformer sampler 所用的标准化方法，确认条件片段化学结构不变；检查输入 token 和参考输出 token 全部存在于本地 prior 的词表中，输入 token 长度不超过 180。仅检查词表/结构，不加载网络执行生成。

当前验证环境：RDKit `{rdBase.rdkitVersion}`，PyTorch `{torch.__version__}`。
{revisions}

| 靶点 | Lib：保留 + 生成重原子 | Link：端 A + 端 B + 生成重原子 | Lib / Link 输入 token |
|---|---:|---:|---:|
{table}

运行方式（从仓库根目录；需已安装 RDKit、PyTorch、Pillow 和 REINVENT4 依赖）：

```bash
/opt/homebrew/Caskroom/miniforge/base/envs/reinvent4/bin/python real-world_dataset/task_design_astra/render_designs.py
```

字体路径 `FONT` 为本机 Arial Unicode；其他机器可改为支持中文的字体。脚本只重写本文件夹内的展示文件，所有候选模型输入字符串仅在内存中用于校验。

## 规范来源

以本地代码为准：

- `REINVENT4/reinvent/runmodes/samplers/libinvent.py`
- `REINVENT4/reinvent/runmodes/samplers/linkinvent.py`
- `REINVENT4/reinvent/chemistry/library_design/attachment_points.py`
- `REINVENT4/reinvent/chemistry/library_design/bond_maker.py`
- `REINVENT4/reinvent/models/transformer/core/vocabulary.py`
- [REINVENT4 官方文档](https://github.com/MolecularAI/REINVENT4/blob/main/contrib/reinvent-doc/index.md)
- [REINVENT4 论文：LibINVENT 骨架与 LinkINVENT 双片段条件](https://pmc.ncbi.nlm.nih.gov/articles/PMC10882833/)

RCSB 数据于 2026-09-22 查阅，CCD 原始下载地址形式为 `https://files.rcsb.org/ligands/download/<CCD>_ideal.sdf`。这些下载文件仅用于核对，没有作为交付输入文件保存。

{sources}
'''
    (OUT/'README.md').write_text(readme,encoding='utf-8')

if __name__=='__main__':
    main()
