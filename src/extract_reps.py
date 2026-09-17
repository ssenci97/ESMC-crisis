#!/usr/bin/env python3
import argparse, datetime, os, resource, time, numpy as np, torch, tqdm
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

MODEL_NAME = "esmc_300m"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
CFG = LogitsConfig(sequence=True, return_embeddings=True)
BATCH_SIZE = 500
START_TIME = time.time()

def parse_fasta(fp):
    r, idx, cid, cseq = [], 0, None, []
    with open(fp) as f:
        for l in f:
            l = l.strip()
            if not l: continue
            if l.startswith(">"):
                if cid: r.append((idx, cid, "".join(cseq))); idx += 1
                cid, cseq = l[1:].split()[0], []
            else: cseq.append(l)
    if cid: r.append((idx, cid, "".join(cseq)))
    return r

def parse_tsv(fp):
    r = []
    with open(fp) as f:
        next(f, None)
        for idx, l in enumerate(f):
            p = l.strip().split("\t")
            if len(p) >= 2: r.append((idx, p[0], p[1]))
    return r

def esmC_extract(seqs, model, device):
    if isinstance(seqs, str): seqs = [seqs]
    tensors = [model.encode(ESMProtein(sequence=s)) for s in seqs]
    return [model.logits(t, CFG).embeddings.squeeze(0)[1:-1].mean(dim=0).cpu().numpy() for t in tensors]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--seqs_tsv")
    g.add_argument("--fasta")
    parser.add_argument("--model_name", default=MODEL_NAME)
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch_size", type=int, default=BATCH_SIZE)
    args = parser.parse_args()
    
    fp = args.fasta if args.fasta else args.seqs_tsv
    records = parse_fasta(args.fasta) if args.fasta else parse_tsv(args.seqs_tsv)
    model = ESMC.from_pretrained(args.model_name).to(DEVICE).eval()
    ids, embs = [], []
    
    with torch.no_grad():
        for i in tqdm.tqdm(range(0, len(records), args.batch_size), desc=f"[PYTHON-INFO] Processing {fp}"):
            b = records[i:i+args.batch_size]
            ids.extend([r[1] for r in b])
            embs.extend(esmC_extract([r[2] for r in b], model, DEVICE))
            
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez(args.out, embeddings=np.stack(embs), ids=np.array(ids))
    print(f"[PYTHON-INFO] date={datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} execution_time_sec={time.time()-START_TIME:.2f} max_memory_mb={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.2f}")
