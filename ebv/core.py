"""Numerical kernels for readout-adaptive event-based velocimetry (RAEBV)."""
import math
import numpy as np
from numba import njit

WEIGHT_EPS = 1e-12
GATE_MIN_CROSS_SPEED_PXMS = 0.5
GATE_MIN_HIST_PROB = 1e-5

@njit(fastmath=True, cache=True)
def build_spatial_lookup(xs, ys, ts, ps, w, h):
    n=len(xs); pix=ys.astype(np.int64)*w+xs.astype(np.int64)
    o=np.argsort(pix)
    xs_s=xs[o]; ys_s=ys[o]; ts_s=ts[o]; ps_s=ps[o]; pix_s=pix[o]
    g=np.zeros((h,w,2),dtype=np.int32)
    if n==0: return xs_s,ys_s,ts_s,ps_s,g
    cp=pix_s[0]; ss=0
    for i in range(1,n):
        if pix_s[i]!=cp:
            gy=cp//w; gx=cp%w
            if 0<=gy<h and 0<=gx<w: g[gy,gx,0]=ss; g[gy,gx,1]=i-ss
            cp=pix_s[i]; ss=i
    gy=cp//w; gx=cp%w
    if 0<=gy<h and 0<=gx<w: g[gy,gx,0]=ss; g[gy,gx,1]=n-ss
    return xs_s,ys_s,ts_s,ps_s,g


@njit(fastmath=True, cache=True)
def pol_curr(cp, pm):
    if pm==0: return cp==1
    if pm==1: return cp==0
    return True


@njit(fastmath=True, cache=True)
def pol_pair(cp, op, pm):
    if pm==0: return op==1
    if pm==1: return op==0
    if pm==2: return op==cp
    return True


def merge_seg(sg):
    if not sg:
        return {'x':np.array([],dtype=np.int32),'y':np.array([],dtype=np.int32),
                't':np.array([],dtype=np.int64),'p':np.array([],dtype=np.int32)}
    return {'x':np.concatenate([s['x'] for s in sg]),
            'y':np.concatenate([s['y'] for s in sg]),
            't':np.concatenate([s['t'] for s in sg]),
            'p':np.concatenate([s['p'] for s in sg])}


def build_triplet(buf, sc, gp):
    bl=len(buf); ws=3*sc+2*gp
    if bl<ws: raise RuntimeError('Not enough segments.')
    s=bl-ws
    return (merge_seg(buf[s:s+sc]),
            merge_seg(buf[s+sc+gp:s+2*sc+gp]),
            merge_seg(buf[s+2*sc+2*gp:s+3*sc+2*gp]))


def _bin_one_layer_for_gated_fusion(ld, w, h, bin_size):
    """Aggregate event velocities with accepted-triplet counts as weights.

    Input: [x, y, u, v, count_meta, triplet_count].
    Output: [x_weighted, y_weighted, u, v, event_count,
             triplet_count_weighted]. Positions and velocities use the same
    weights. The last column is metadata, not a statistical uncertainty.
    """
    if ld is None or len(ld) == 0:
        return np.zeros((0, 6), dtype=np.float64)

    x = ld[:, 0]; y = ld[:, 1]
    u = ld[:, 2]; v = ld[:, 3]
    tc = ld[:, 5]

    valid = (np.isfinite(x) & np.isfinite(y) &
             np.isfinite(u) & np.isfinite(v) &
             np.isfinite(tc))
    if not np.any(valid):
        return np.zeros((0, 6), dtype=np.float64)

    x = x[valid]; y = y[valid]
    u = u[valid]; v = v[valid]
    tc = tc[valid]

    grid_w = (w + bin_size - 1) // bin_size
    grid_h = (h + bin_size - 1) // bin_size

    gx = (x // bin_size).astype(np.int32)
    gy = (y // bin_size).astype(np.int32)
    inside = (gx >= 0) & (gx < grid_w) & (gy >= 0) & (gy < grid_h)
    if not np.any(inside):
        return np.zeros((0, 6), dtype=np.float64)

    gx = gx[inside]; gy = gy[inside]
    x = x[inside]; y = y[inside]
    u = u[inside]; v = v[inside]
    tc = tc[inside]

    wt = np.maximum(tc, 1.0)

    n   = np.zeros((grid_h, grid_w), dtype=np.float64)  # raw vector count
    sw  = np.zeros((grid_h, grid_w), dtype=np.float64)  # triplet-count weight sum
    sx  = np.zeros((grid_h, grid_w), dtype=np.float64)
    sy  = np.zeros((grid_h, grid_w), dtype=np.float64)
    su  = np.zeros((grid_h, grid_w), dtype=np.float64)
    sv  = np.zeros((grid_h, grid_w), dtype=np.float64)
    stc = np.zeros((grid_h, grid_w), dtype=np.float64)

    np.add.at(n,   (gy, gx), 1.0)
    np.add.at(sw,  (gy, gx), wt)
    np.add.at(sx,  (gy, gx), wt * x)
    np.add.at(sy,  (gy, gx), wt * y)
    np.add.at(su,  (gy, gx), wt * u)
    np.add.at(sv,  (gy, gx), wt * v)
    np.add.at(stc, (gy, gx), wt * tc)

    rows, cols = np.where(sw > 0.0)
    if len(rows) == 0:
        return np.zeros((0, 6), dtype=np.float64)

    ww = sw[rows, cols]
    nn = n[rows, cols]

    out = np.column_stack((
        sx[rows, cols] / ww,
        sy[rows, cols] / ww,
        su[rows, cols] / ww,
        sv[rows, cols] / ww,
        nn,
        stc[rows, cols] / ww,
    ))
    ok = np.isfinite(out).all(axis=1) & (out[:, 4] > 0.0)
    return out[ok]


def _binned_layer_to_maps(binned, w, h, bin_size):
    """Convert a layer-wise binned vector table into grid maps."""
    grid_w = (w + bin_size - 1) // bin_size
    grid_h = (h + bin_size - 1) // bin_size
    maps = [np.zeros((grid_h, grid_w), dtype=np.float64) for _ in range(6)]
    mask = np.zeros((grid_h, grid_w), dtype=bool)

    if binned is None or len(binned) == 0:
        return maps, mask

    gx = (binned[:, 0] // bin_size).astype(np.int32)
    gy = (binned[:, 1] // bin_size).astype(np.int32)
    valid = (gx >= 0) & (gx < grid_w) & (gy >= 0) & (gy < grid_h) & np.isfinite(binned).all(axis=1)
    gx = gx[valid]; gy = gy[valid]
    bb = binned[valid]

    if len(bb) == 0:
        return maps, mask

    for j in range(6):
        maps[j][gy, gx] = bb[:, j]
    mask[gy, gx] = bb[:, 4] > WEIGHT_EPS
    return maps, mask


def fuse_layers_speed_gated_bin(lr, w, h, bin_size, speed_gate_pxs):
    """Bin each layer and select one velocity per cell.

    Cell velocities are weighted by accepted-triplet counts. Where both
    layers are valid, select layer 1 if the layer-0 speed is below the gate,
    otherwise select layer 0. Where only one layer is valid, retain it.

    Output columns: [x, y, u, v, layer_id].
    """
    if lr is None or len(lr) == 0:
        return None

    if len(lr) == 1:
        l0 = _bin_one_layer_for_gated_fusion(lr[0], w, h, bin_size)
        if l0 is None or len(l0) == 0:
            return None
        return np.column_stack((l0[:, 0], l0[:, 1], l0[:, 2], l0[:, 3],
                                np.zeros(len(l0), dtype=np.float64)))

    if len(lr) != 2:
        raise NotImplementedError(
            "speed-gated fusion currently supports exactly 1 or 2 layers. "
            "N>2 is intentionally not implemented."
        )

    l0 = _bin_one_layer_for_gated_fusion(lr[0], w, h, bin_size)
    l1 = _bin_one_layer_for_gated_fusion(lr[1], w, h, bin_size)

    m0, mask0 = _binned_layer_to_maps(l0, w, h, bin_size)
    m1, mask1 = _binned_layer_to_maps(l1, w, h, bin_size)

    union = mask0 | mask1
    rows, cols = np.where(union)
    n = len(rows)
    if n == 0:
        return None

    out = np.full((n, 5), np.nan, dtype=np.float64)

    has0 = mask0[rows, cols]
    has1 = mask1[rows, cols]

    only0 = has0 & (~has1)
    only1 = has1 & (~has0)
    both  = has0 & has1

    layer_label = np.full(mask0.shape, -1, dtype=np.int8)
    layer_label[mask0 & (~mask1)] = 0
    layer_label[mask1 & (~mask0)] = 1
    both_rr = both_cc = both_use_layer1 = None
    if np.any(both):
        both_rr = rows[both]
        both_cc = cols[both]
        both_u0 = m0[2][both_rr, both_cc]
        both_v0 = m0[3][both_rr, both_cc]
        both_s0 = np.sqrt(both_u0*both_u0 + both_v0*both_v0)
        both_use_layer1 = both_s0 < speed_gate_pxs
        layer_label[both_rr, both_cc] = np.where(both_use_layer1, 1, 0).astype(np.int8)

    if np.any(only0):
        rr = rows[only0]
        cc = cols[only0]
        out[only0, 0] = m0[0][rr, cc]
        out[only0, 1] = m0[1][rr, cc]
        out[only0, 2] = m0[2][rr, cc]
        out[only0, 3] = m0[3][rr, cc]
        out[only0, 4] = 0.0

    if np.any(only1):
        rr = rows[only1]
        cc = cols[only1]
        out[only1, 0] = m1[0][rr, cc]
        out[only1, 1] = m1[1][rr, cc]
        out[only1, 2] = m1[2][rr, cc]
        out[only1, 3] = m1[3][rr, cc]
        out[only1, 4] = 1.0

    if np.any(both):
        rr = both_rr
        cc = both_cc
        use_layer1 = both_use_layer1
        both_idx = np.where(both)[0]

        if np.any(use_layer1):
            idx = both_idx[use_layer1]
            rr_l = rr[use_layer1]
            cc_l = cc[use_layer1]
            out[idx, 0] = m1[0][rr_l, cc_l]
            out[idx, 1] = m1[1][rr_l, cc_l]
            out[idx, 2] = m1[2][rr_l, cc_l]
            out[idx, 3] = m1[3][rr_l, cc_l]
            out[idx, 4] = 1.0

        if np.any(~use_layer1):
            idx = both_idx[~use_layer1]
            rr_h = rr[~use_layer1]
            cc_h = cc[~use_layer1]
            out[idx, 0] = m0[0][rr_h, cc_h]
            out[idx, 1] = m0[1][rr_h, cc_h]
            out[idx, 2] = m0[2][rr_h, cc_h]
            out[idx, 3] = m0[3][rr_h, cc_h]
            out[idx, 4] = 0.0

    ok = np.isfinite(out).all(axis=1)
    if not np.any(ok):
        return None
    return out[ok]


def _smooth_hist_prob(p, window):
    p = np.asarray(p, dtype=np.float64)
    if window is None or window <= 1:
        return p.copy()
    w = int(window)
    if w < 1:
        return p.copy()
    kernel = np.ones(w, dtype=np.float64) / float(w)
    ps = np.convolve(p, kernel, mode="same")
    s = ps.sum()
    if s > 0:
        ps /= s
    return ps


def _find_hist_intersection(speed_centers, p0, p1):
    """Find the meaningful low-to-high layer speed crossover.

    Returns:
        center_pxms, method_string

    The selected crossing is the first p1-dominant -> p0-dominant crossing
    after the layer-1 peak and above GATE_MIN_CROSS_SPEED_PXMS. If no such
    crossing exists, use the minimum absolute difference after that search start
    as a conservative fallback and mark it in the method string.
    """
    c = np.asarray(speed_centers, dtype=np.float64)
    p0 = np.asarray(p0, dtype=np.float64)
    p1 = np.asarray(p1, dtype=np.float64)
    if len(c) < 2 or len(p0) != len(c) or len(p1) != len(c):
        raise RuntimeError("Invalid histogram arrays for gate intersection search.")
    if np.sum(p0) <= 0 or np.sum(p1) <= 0:
        raise RuntimeError("Empty layer histogram; cannot calibrate speed gate.")

    diff = p0 - p1
    min_idx = int(np.searchsorted(c, GATE_MIN_CROSS_SPEED_PXMS))
    peak1_idx = int(np.argmax(p1))
    start = max(min_idx, peak1_idx + 1, 1)
    start = min(start, len(c) - 1)

    # Preferred crossing: p1>p0 below, p0>=p1 above.
    for i in range(start, len(c)):
        if max(p0[i-1], p1[i-1], p0[i], p1[i]) < GATE_MIN_HIST_PROB:
            continue
        if diff[i-1] <= 0.0 and diff[i] >= 0.0:
            d0 = diff[i-1]
            d1 = diff[i]
            x0 = c[i-1]
            x1 = c[i]
            if abs(d1 - d0) > 1e-15:
                x = x0 + (0.0 - d0) * (x1 - x0) / (d1 - d0)
            else:
                x = 0.5 * (x0 + x1)
            return float(x), "p1_to_p0_crossing"

    # Fallback: closest p0/p1 probability after start where the histograms are non-negligible.
    valid = (c >= GATE_MIN_CROSS_SPEED_PXMS) & ((p0 + p1) >= GATE_MIN_HIST_PROB)
    valid[:start] = False
    if not np.any(valid):
        valid = (c >= GATE_MIN_CROSS_SPEED_PXMS)
    idxs = np.where(valid)[0]
    j = idxs[np.argmin(np.abs(diff[idxs]))]
    return float(c[j]), "fallback_min_abs_hist_difference"


@njit(fastmath=True, cache=True)
def fit_vel_3pt(t0, z0, t1, z1, t2, z2):
    """Returns (slope [px/s], RSS_z [px^2], sum_dt2 [us^2])"""
    tm=(t0+t1+t2)/3.0; zm=(z0+z1+z2)/3.0
    stt=(t0-tm)**2+(t1-tm)**2+(t2-tm)**2
    stz=(t0-tm)*(z0-zm)+(t1-tm)*(z1-zm)+(t2-tm)*(z2-zm)
    slope=stz/(stt+1e-12)
    intercept=zm-slope*tm
    e0=z0-(intercept+slope*t0)
    e1=z1-(intercept+slope*t1)
    e2=z2-(intercept+slope*t2)
    rss=e0*e0+e1*e1+e2*e2
    return slope*1e6, rss, stt


@njit(fastmath=True, cache=True)
def process_triplet_v8(
    curr_xs, curr_ys, curr_ts, curr_ps,
    prev_xs, prev_ys, prev_ts, prev_ps, prev_grid,
    pprev_xs, pprev_ys, pprev_ts, pprev_ps, pprev_grid,
    width, height,
    x_min_offset, x_max_offset, y_min_offset, y_max_offset,
    polarity_mode,
    min_mag, max_mag,
    pprev_t_min, pprev_t_max,
    trajectory_tol_px,
):
    """
    Match collinear event triplets and average their velocities.

    For each current event:
      1. For each prev event in search radius, define a constant-velocity trajectory
      2. Gather pprev candidates from a coarse bounding box around that trajectory
      3. Accept a pprev event only if its radial residual from the trajectory
         at its own timestamp is <= trajectory_tol_px
      4. Compute velocity (3-point linear fit) and accumulate uniformly

    trajectory_tol_px specifies the radial trajectory tolerance in pixels.

    Output columns: [x, y, vx, vy, count, count]
    """
    nc = len(curr_xs)
    kept = np.zeros((nc, 6), dtype=np.float64)
    kc = 0

    # Coarse integer padding for candidate lookup. Exact acceptance below uses
    # radial distance to the trajectory, not this rectangular bounding box.
    pad = int(math.ceil(trajectory_tol_px))
    if pad < 1:
        pad = 1
    traj_tol2 = trajectory_tol_px * trajectory_tol_px

    for i in range(nc):
        cx = float(curr_xs[i]); cy = float(curr_ys[i])
        ct = float(curr_ts[i]); cp = curr_ps[i]
        if not pol_curr(cp, polarity_mode):
            continue

        # Accumulators for uniform average
        sum_vx = 0.0
        sum_vy = 0.0
        cnt    = 0

        xmn = max(0, int(cx + x_min_offset))
        xmx = min(width - 1, int(cx + x_max_offset))
        ymn = max(0, int(cy + y_min_offset))
        ymx = min(height - 1, int(cy + y_max_offset))

        for ax in range(xmn, xmx + 1):
            for ay in range(ymn, ymx + 1):
                si = int(prev_grid[ay, ax, 0])
                cn = int(prev_grid[ay, ax, 1])
                if cn == 0:
                    continue
                for k in range(si, si + cn):
                    pt = float(prev_ts[k]); pp = prev_ps[k]
                    if not pol_pair(cp, pp, polarity_mode):
                        continue
                    d1 = ct - pt
                    if d1 <= 0.0:
                        continue
                    dx_ = cx - float(prev_xs[k])
                    dy_ = cy - float(prev_ys[k])
                    cvx = (dx_ / d1) * 1e6
                    cvy = (dy_ / d1) * 1e6
                    fm = math.sqrt(cvx * cvx + cvy * cvy)
                    if not (min_mag <= fm <= max_mag):
                        continue

                    # Coarse pprev candidate lookup: bounding box of the
                    # curr-prev trajectory over the pprev time range.
                    # Exact radial trajectory consistency is checked per event below.
                    pxa = cx - cvx * ((ct - float(pprev_t_min)) * 1e-6)
                    pxb = cx - cvx * ((ct - float(pprev_t_max)) * 1e-6)
                    pya = cy - cvy * ((ct - float(pprev_t_min)) * 1e-6)
                    pyb = cy - cvy * ((ct - float(pprev_t_max)) * 1e-6)
                    sxn = max(0, int(math.floor(min(pxa, pxb) - pad)))
                    sxx = min(width - 1, int(math.ceil(max(pxa, pxb) + pad)))
                    syn = max(0, int(math.floor(min(pya, pyb) - pad)))
                    syx = min(height - 1, int(math.ceil(max(pya, pyb) + pad)))
                    if sxn > sxx or syn > syx:
                        continue

                    for xj in range(sxn, sxx + 1):
                        for yj in range(syn, syx + 1):
                            sj = int(pprev_grid[yj, xj, 0])
                            cj = int(pprev_grid[yj, xj, 1])
                            if cj == 0:
                                continue
                            for m in range(sj, sj + cj):
                                ppt = float(pprev_ts[m])
                                ppp = pprev_ps[m]
                                if not pol_pair(cp, ppp, polarity_mode):
                                    continue
                                d2 = pt - ppt
                                if d2 <= 0.0:
                                    continue
                                db = ct - ppt
                                if db <= 0.0:
                                    continue

                                ppx = float(pprev_xs[m]); ppy = float(pprev_ys[m])

                                # Fine strict check: the pprev event must lie within
                                # the radial tolerance from the curr-prev constant-
                                # velocity trajectory at its own timestamp.
                                x_pred = cx - cvx * ((ct - ppt) * 1e-6)
                                y_pred = cy - cvy * ((ct - ppt) * 1e-6)
                                rx = ppx - x_pred
                                ry = ppy - y_pred
                                if rx * rx + ry * ry > traj_tol2:
                                    continue

                                # Velocity estimate — 3-point linear fit
                                # (uses all three events; variance reduction
                                #  via non-uniform timing)
                                pvx = float(prev_xs[k]); pvy = float(prev_ys[k])
                                vxt, _, _ = fit_vel_3pt(ppt, ppx, pt, pvx, ct, cx)
                                vyt, _, _ = fit_vel_3pt(ppt, ppy, pt, pvy, ct, cy)

                                if (not np.isfinite(vxt)) or (not np.isfinite(vyt)):
                                    continue
                                fm2 = math.sqrt(vxt * vxt + vyt * vyt)
                                if not (min_mag <= fm2 <= max_mag):
                                    continue

                                # Uniform accumulation — no score weighting
                                sum_vx += vxt
                                sum_vy += vyt
                                cnt    += 1

        # Mean velocity over accepted triplets.
        if cnt == 0:
            continue

        vx = sum_vx / float(cnt)
        vy = sum_vy / float(cnt)
        if (not np.isfinite(vx)) or (not np.isfinite(vy)):
            continue
        fm = math.sqrt(vx * vx + vy * vy)
        if not (min_mag <= fm <= max_mag):
            continue

        kept[kc, 0] = cx
        kept[kc, 1] = cy
        kept[kc, 2] = vx
        kept[kc, 3] = vy
        kept[kc, 4] = float(cnt)   # count_meta = triplet count before binning
        kept[kc, 5] = float(cnt)
        kc += 1

    return kept[:kc]
