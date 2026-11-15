"""CRNN segment classifier for long-billed hermit songs.

Whole-song spectrogram -> 2D CNN -> BiGRU over time -> per-frame logits,
averaged within each annotated segment window -> 6-class label per segment.
Grouped K-fold CV by consensus song type (lek.song.type.CONS).
"""
import sys, json
import numpy as np, pandas as pd, torch, torch.nn as nn
from scipy.io import wavfile
from scipy.signal import resample_poly
from sklearn.model_selection import GroupKFold
from sklearn.metrics import cohen_kappa_score, confusion_matrix, accuracy_score

S = sys.argv[1] if len(sys.argv) > 1 else "."
GROUP = sys.argv[2] if len(sys.argv) > 2 else "lek.song.type.CONS"
torch.manual_seed(0); np.random.seed(0)
dev = "cuda"
SR, HOP, NFFT = 44100, 22, 256          # ~0.5 ms hop -> ~14 frames per 7 ms segment
FMIN, FMAX = 1000, 12000
CLASSES = ["d", "f", "m", "p", "s", "u"]

seg = pd.read_csv(f"{S}/segments.csv")
songs = sorted(seg["sound.files"].unique())
freqs = np.linspace(0, SR / 2, NFFT // 2 + 1); fsel = (freqs >= FMIN) & (freqs <= FMAX)

data = []
for s in songs:
    sr, wav = wavfile.read(f"{S}/wav/{s}.wav"); wav = wav.astype(np.float64)
    if sr != SR: wav = resample_poly(wav, SR, sr)
    wav = torch.tensor(wav - wav.mean()); wav = wav / (wav.abs().max() + 1e-9)
    st = torch.stft(wav, NFFT, HOP, window=torch.hann_window(NFFT, dtype=wav.dtype), return_complex=True).abs() ** 2
    sp = torch.log10(st[torch.tensor(fsel)] + 1e-10)
    sp = (sp - sp.mean()) / (sp.std() + 1e-9)
    d = seg[seg["sound.files"] == s].sort_values("pos")
    nfr = sp.shape[1]
    win = [(min(int(a * SR / HOP), nfr - 1), max(min(int(np.ceil(b * SR / HOP)), nfr), int(a * SR / HOP) + 1)) for a, b in zip(d.start, d.end)]
    data.append(dict(song=s, spec=sp.float(), win=win,
                     y=torch.tensor([CLASSES.index(l) for l in d.label6]),
                     group=d[GROUP].iloc[0], lek=d.lek.iloc[0]))
NF = data[0]["spec"].shape[0]
print(f"{len(data)} songs, {NF} freq bins", flush=True)


class CRNN(nn.Module):
    def __init__(s, nf, nc=6, h=96):
        super().__init__()
        def blk(i, o): return nn.Sequential(nn.Conv2d(i, o, 3, padding=1), nn.BatchNorm2d(o), nn.GELU(),
                                            nn.Conv2d(o, o, 3, padding=1), nn.BatchNorm2d(o), nn.GELU(),
                                            nn.MaxPool2d((2, 1)))
        s.cnn = nn.Sequential(blk(1, 24), blk(24, 48), blk(48, 64), nn.Dropout2d(0.2))
        s.proj = nn.Linear(64 * (nf // 8), h)
        s.rnn = nn.GRU(h, h, num_layers=2, batch_first=True, bidirectional=True, dropout=0.3)
        s.out = nn.Sequential(nn.Dropout(0.3), nn.Linear(2 * h, nc))

    def forward(s, x):                         # x: B,F,T
        z = s.cnn(x.unsqueeze(1))              # B,C,F',T
        B, C, F, T = z.shape
        z = s.proj(z.permute(0, 3, 1, 2).reshape(B, T, C * F))
        z, _ = s.rnn(z)
        return s.out(z)                        # B,T,nc


def seg_logits(fr, wins):                      # fr: T,nc  -> nseg,nc
    return torch.stack([fr[a:b].mean(0) for a, b in wins])


def augment(sp):
    sp = sp.clone()
    sp += 0.15 * torch.randn_like(sp)                       # noise
    f0 = np.random.randint(0, sp.shape[0] - 10); sp[f0:f0 + np.random.randint(0, 10)] = 0  # freq mask
    sh = np.random.randint(-4, 5); sp = torch.roll(sp, sh, 0)  # small pitch shift
    return sp * np.random.uniform(0.8, 1.2)


def run_batch(model, items, train):
    T = max(it["spec"].shape[1] for it in items)
    X = torch.zeros(len(items), NF, T)
    for i, it in enumerate(items):
        sp = augment(it["spec"]) if train else it["spec"]
        X[i, :, :sp.shape[1]] = sp
    fr = model(X.to(dev))
    return torch.cat([seg_logits(fr[i], it["win"]) for i, it in enumerate(items)]), torch.cat([it["y"] for it in items]).to(dev)


groups = np.array([d["group"] for d in data])
counts = np.bincount(torch.cat([d["y"] for d in data]).numpy(), minlength=6)
cw = torch.tensor((counts.sum() / counts / 6) ** 0.5, dtype=torch.float32, device=dev)
preds = {}
EPOCHS, BS = 60, 16
for k, (tr, te) in enumerate(GroupKFold(n_splits=5).split(np.zeros(len(data)), groups=groups)):
    model = CRNN(NF).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-3)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=2e-3, total_steps=EPOCHS * int(np.ceil(len(tr) / BS)))
    lossf = nn.CrossEntropyLoss(weight=cw, label_smoothing=0.05)
    for ep in range(EPOCHS):
        model.train(); perm = np.random.permutation(tr)
        for i in range(0, len(perm), BS):
            lg, y = run_batch(model, [data[j] for j in perm[i:i + BS]], True)
            loss = lossf(lg, y); opt.zero_grad(); loss.backward(); opt.step(); sched.step()
    model.eval()
    with torch.no_grad():
        for j in te:
            lg, _ = run_batch(model, [data[j]], False)
            preds[data[j]["song"]] = torch.softmax(lg, -1).cpu().numpy()
    yt = torch.cat([data[j]["y"] for j in te]).numpy(); yp = np.concatenate([preds[data[j]["song"]].argmax(1) for j in te])
    print(f"fold {k}: acc={accuracy_score(yt, yp):.3f} kappa={cohen_kappa_score(yt, yp):.3f}", flush=True)

yt = torch.cat([d["y"] for d in data]).numpy(); yp = np.concatenate([preds[d["song"]].argmax(1) for d in data])
print(f"\nOVERALL ({GROUP} held out): acc={accuracy_score(yt, yp):.3f} kappa={cohen_kappa_score(yt, yp):.3f}")
print("rows=true cols=pred", CLASSES); print(confusion_matrix(yt, yp))
out = []
for d in data:
    for i, p in enumerate(preds[d["song"]]):
        out.append(dict(song=d["song"], pos=i + 1, true=CLASSES[d["y"][i]], pred=CLASSES[p.argmax()], conf=p.max(), **{f"p_{c}": p[ci] for ci, c in enumerate(CLASSES)}))
pd.DataFrame(out).to_csv(f"{S}/cv_predictions_{GROUP}.csv", index=False)
