"""Direct cell assignment and equal-window population moments."""
from pathlib import Path
import numpy as np
import pandas as pd


def grid_shape(c):
    step = c["cell_size_px"]
    return tuple((c["roi"][k]+step-1)//step for k in ("height", "width"))


def direct_grid(vec, c):
    shape = grid_shape(c)
    out = np.full((*shape, 2), np.nan)
    if not len(vec):
        return out
    x = (vec[:, 0]//c["cell_size_px"]).astype(int)
    y = (vec[:, 1]//c["cell_size_px"]).astype(int)
    ids = y*shape[1]+x
    if len(np.unique(ids)) != len(vec) or not np.isfinite(vec[:, :4]).all():
        raise ValueError("Expected one finite MAIN vector per cell")
    if not np.all((x >= 0) & (x < shape[1]) & (y >= 0) & (y < shape[0])):
        raise ValueError("Vector outside the configured cell grid")
    out[y, x] = vec[:, 2:4]
    return out


class Moments:
    def __init__(self, shape):
        self.count = np.zeros(shape, dtype=np.int64)
        self.mean = np.zeros((*shape, 2))
        self.m2 = np.zeros((*shape, 2))
        self.cross = np.zeros(shape)

    def add(self, a):
        ok = np.isfinite(a).all(axis=2)
        n = self.count[ok]+1
        d = a[ok]-self.mean[ok]
        self.count[ok] = n
        self.mean[ok] += d/n[:, None]
        e = a[ok]-self.mean[ok]
        self.m2[ok] += d*e
        self.cross[ok] += d[:, 0]*e[:, 1]

    def merge(self, src):
        na, nb = self.count.astype(float), src.count.astype(float)
        total = na+nb
        valid = total > 0
        ratio = np.zeros_like(total)
        cross_weight = np.zeros_like(total)
        ratio[valid] = nb[valid]/total[valid]
        cross_weight[valid] = na[valid]*nb[valid]/total[valid]
        delta = src.mean-self.mean
        self.mean += delta*ratio[:, :, None]
        self.m2 += src.m2+delta*delta*cross_weight[:, :, None]
        self.cross += src.cross+delta[:, :, 0]*delta[:, :, 1]*cross_weight
        self.count += src.count

    def table(self, c, si=False):
        step = c["cell_size_px"]
        yy, xx = np.indices(self.count.shape)
        n = np.maximum(self.count, 1)
        rms = np.sqrt(np.maximum(0.0, self.m2/n[:, :, None]))
        fac = c["spatial_scale_m_per_px"] if si else 1.0
        d = pd.DataFrame(dict(
            x=((xx+.5)*step).ravel()*fac, y=((yy+.5)*step).ravel()*fac,
            u=self.mean[:, :, 0].ravel()*fac, v=self.mean[:, :, 1].ravel()*fac,
            u_rms=rms[:, :, 0].ravel()*fac, v_rms=rms[:, :, 1].ravel()*fac,
            uv_stress=(-self.cross/n).ravel()*fac**2,
            ke=(.5*(rms**2).sum(axis=2)).ravel()*fac**2,
            count=self.count.ravel()))
        return d[d["count"] >= c["min_valid_samples"]].copy()

    def save(self, path):
        path = Path(path)
        temp = path.with_suffix(".tmp")
        with temp.open("wb") as fh:
            np.savez_compressed(fh, count=self.count, mean=self.mean, m2=self.m2, cross=self.cross)
        temp.replace(path)

    def restore(self, path):
        with np.load(path, allow_pickle=False) as data:
            for key in ("count", "mean", "m2", "cross"):
                a = data[key]
                if a.shape != getattr(self, key).shape:
                    raise ValueError(f"Wrong moment shape: {key}")
                setattr(self, key, a.copy())
