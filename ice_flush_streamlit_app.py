
# -*- coding: utf-8 -*-
import re, io, zipfile, math
import pandas as pd
import numpy as np
import streamlit as st
import plotly.graph_objects as go

st.set_page_config(page_title="冰浆冲洗可视化核算", layout="wide")

def read_csv_bytes(b):
    for enc in ("utf-8-sig","utf-8","gbk","gb2312"):
        try:
            s=b.decode(enc)
            return pd.read_csv(io.StringIO(s))
        except Exception:
            pass
    return pd.read_csv(io.BytesIO(b), encoding="gbk", errors="ignore")

def read_table_bytes(b, name):
    if name.lower().endswith(".csv"):
        return read_csv_bytes(b)
    return pd.read_excel(io.BytesIO(b))

def infer_date_from_names(names):
    for n in names:
        m=re.search(r"(20\d{6})", n)
        if m:
            s=m.group(1)
            return f"{s[:4]}-{s[4:6]}-{s[6:8]}"
    return None

def pick_by_keywords(names, keywords):
    # 优先：同时含所有关键字
    best=[]
    for n in names:
        if all(k in n for k in keywords):
            best.append(n)
    if best:
        return best[0]
    # 退化：含任一关键字
    for n in names:
        if any(k in n for k in keywords):
            return n
    return None

def pick_numeric_col(df, prefer_contains=None):
    cols=[c for c in df.columns if c!="TIME"]
    if prefer_contains:
        for key in prefer_contains:
            for c in cols:
                if key in str(c):
                    return c
    for c in cols:
        try:
            pd.to_numeric(df[c].iloc[:20], errors="raise")
            return c
        except Exception:
            continue
    return cols[0] if cols else None

def parse_dt(date_str, time_series):
    t=time_series.astype(str).str.strip()
    return pd.to_datetime(date_str + " " + t, errors="coerce")

def merge_site(date_str, sludge, flow, turb, ec, temp, tol_sec=8):
    def prep(df):
        d=df.copy()
        if "TIME" not in d.columns:
            d = d.rename(columns={d.columns[0]:"TIME"})
        d["dt"]=parse_dt(date_str, d["TIME"])
        d=d.dropna(subset=["dt"]).sort_values("dt")
        return d

    s=prep(sludge); f=prep(flow); t=prep(turb); e=prep(ec); p=prep(temp)
    s_col=pick_numeric_col(s, ["测值","值","浓度"])
    f_col=pick_numeric_col(f, ["瞬时流量","流量"])
    t_col=pick_numeric_col(t, ["浊度"])
    e_col=pick_numeric_col(e, ["测量值","电导"])
    p_col=pick_numeric_col(p, ["温度"])

    s=s[["dt",s_col]].rename(columns={s_col:"sludge_gL"})
    f=f[["dt",f_col]].rename(columns={f_col:"Q_m3h"})
    t=t[["dt",t_col]].rename(columns={t_col:"turb_NTU"})
    e=e[["dt",e_col]].rename(columns={e_col:"EC_mScm"})
    p=p[["dt",p_col]].rename(columns={p_col:"temp_C"})
    for d,c in [(s,"sludge_gL"),(f,"Q_m3h"),(t,"turb_NTU"),(e,"EC_mScm"),(p,"temp_C")]:
        d[c]=pd.to_numeric(d[c], errors="coerce")

    tol=pd.Timedelta(seconds=tol_sec)
    out=pd.merge_asof(s,f,on="dt",direction="nearest",tolerance=tol)
    out=pd.merge_asof(out,t,on="dt",direction="nearest",tolerance=tol)
    out=pd.merge_asof(out,e,on="dt",direction="nearest",tolerance=tol)
    out=pd.merge_asof(out,p,on="dt",direction="nearest",tolerance=tol)
    return out

def compute_c0(df, ec_thresh=0.6, turb_thresh=5):
    filt=(df["EC_mScm"]<=ec_thresh)&(df["turb_NTU"]<=turb_thresh)
    if filt.sum()>=10:
        return float(df.loc[filt,"sludge_gL"].median())
    return float(df["sludge_gL"].median())

def integrate(df, c0, ec_ice_thresh=5.0, ice_vol=None, total_vol=None, expected_sec=10):
    d=df.copy()
    d["dt_s"]=d["dt"].diff().dt.total_seconds().fillna(0)
    d["dt_h"]=d["dt_s"]/3600
    d["dV"]=d["Q_m3h"].abs()*d["dt_h"]
    d["C_ex"]=(d["sludge_gL"]-c0).clip(lower=0)
    d["dM"]=d["C_ex"]*d["dV"]
    d["cumV"]=d["dV"].cumsum()
    d["cumM"]=d["dM"].cumsum()
    # data quality
    duration=max((d["dt"].iloc[-1]-d["dt"].iloc[0]).total_seconds(),1)
    missing_ratio=float(np.sum(np.clip(d["dt_s"].values[1:]-expected_sec,0,None))/duration)
    q0_ratio=float((d["Q_m3h"]==0).mean())

    ice_start_idx = d.index[d["EC_mScm"]>=ec_ice_thresh]
    ice_start_idx = int(ice_start_idx[0]) if len(ice_start_idx) else None
    ice_end_idx=None
    total_end_idx=None

    if total_vol is not None:
        idxs=d.index[d["cumV"]>=total_vol]
        total_end_idx=int(idxs[0]) if len(idxs) else None

    cumM_before=0.0
    if ice_start_idx is not None and ice_vol is not None:
        cumV_at=float(d.loc[ice_start_idx,"cumV"])
        cumM_before=float(d.loc[ice_start_idx-1,"cumM"]) if ice_start_idx>0 else 0.0
        d["cumV_ice"]=d["cumV"]-cumV_at+d["dV"]
        idxs=d.index[(d.index>=ice_start_idx)&(d["cumV_ice"]>=ice_vol)]
        ice_end_idx=int(idxs[0]) if len(idxs) else None

    M_total=float(d.loc[total_end_idx,"cumM"]) if total_end_idx is not None else float(d["cumM"].iloc[-1])
    M_ice=float(d.loc[ice_end_idx,"cumM"]-cumM_before) if ice_end_idx is not None else float("nan")

    return d, dict(
        c0=c0,
        ice_start_time=d.loc[ice_start_idx,"dt"] if ice_start_idx is not None else None,
        ice_end_time=d.loc[ice_end_idx,"dt"] if ice_end_idx is not None else None,
        total_end_time=d.loc[total_end_idx,"dt"] if total_end_idx is not None else None,
        M_total_kg=M_total,
        M_ice_kg=M_ice,
        q0_ratio=q0_ratio,
        missing_ratio=missing_ratio
    )

def plot_timeseries(d):
    fig=go.Figure()
    fig.add_trace(go.Scatter(x=d["dt"], y=d["Q_m3h"], name="Q(m³/h)", yaxis="y1"))
    fig.add_trace(go.Scatter(x=d["dt"], y=d["sludge_gL"], name="污泥C(g/L)", yaxis="y2"))
    fig.add_trace(go.Scatter(x=d["dt"], y=d["turb_NTU"], name="浊度(NTU)", yaxis="y3"))
    fig.add_trace(go.Scatter(x=d["dt"], y=d["EC_mScm"], name="电导率(mS/cm)", yaxis="y4"))
    fig.update_layout(
        height=520,
        xaxis=dict(domain=[0.06,0.98]),
        yaxis=dict(title="Q", side="left", position=0.02),
        yaxis2=dict(title="C", overlaying="y", side="right", position=0.98),
        yaxis3=dict(title="Turb", overlaying="y", side="right", position=0.94, showgrid=False),
        yaxis4=dict(title="EC", overlaying="y", side="right", position=0.90, showgrid=False),
        legend=dict(orientation="h", y=1.02, x=0),
        margin=dict(l=40,r=20,t=30,b=40),
    )
    return fig

def plot_cumulative(d, res):
    fig=go.Figure()
    fig.add_trace(go.Scatter(x=d["dt"], y=d["cumV"], name="CumV(m³)", yaxis="y1"))
    fig.add_trace(go.Scatter(x=d["dt"], y=d["cumM"], name="CumM(kg)", yaxis="y2"))
    for t in [res["ice_start_time"], res["ice_end_time"], res["total_end_time"]]:
        if t is None: 
            continue
        fig.add_vline(x=t, line_width=2, line_dash="dot")
    fig.update_layout(
        height=520,
        xaxis=dict(domain=[0.06,0.98]),
        yaxis=dict(title="CumV", side="left", position=0.02),
        yaxis2=dict(title="CumM", overlaying="y", side="right", position=0.98),
        legend=dict(orientation="h", y=1.02, x=0),
        margin=dict(l=40,r=20,t=30,b=40),
    )
    return fig

st.title("冰浆冲洗｜去除量核算（Python前端网页）")
st.caption("上传一个点位的 zip（里面含 污泥/流量/浊度/电导/温度 文件），设置参数后自动计算 M_total / M_ice 并可视化。")

uploaded = st.file_uploader("上传点位 zip（或包含该点位 CSV/XLSX 的 zip）", type=["zip"])
colA, colB, colC, colD = st.columns(4)
with colA:
    date_str = st.text_input("日期(YYYY-MM-DD)", value="")
with colB:
    tol_sec = st.number_input("时间对齐容差(s)", min_value=1, max_value=60, value=8)
with colC:
    ec_ice = st.number_input("电导阈值(冰浆开始, mS/cm)", min_value=0.0, value=5.0, step=0.1)
with colD:
    expected_sec = st.number_input("期望采样间隔(s)", min_value=1, max_value=120, value=10)

col1, col2, col3 = st.columns(3)
with col1:
    ice_vol = st.number_input("冰浆耗水(m³)用于切段", min_value=0.0, value=0.0, step=0.1)
with col2:
    total_vol = st.number_input("总耗水(m³)用于切段(可为0表示用数据末端)", min_value=0.0, value=0.0, step=0.1)
with col3:
    report_mass = st.number_input("报告质量(kg,用于校准k，可为0)", min_value=0.0, value=0.0, step=0.1)

if uploaded is None:
    st.info("先上传一个 zip 文件。")
    st.stop()

# read zip
with zipfile.ZipFile(uploaded, "r") as z:
    names = z.namelist()

auto_date = infer_date_from_names(names)
if not date_str:
    date_str = auto_date or "2025-01-01"
    st.write(f"自动识别日期：{date_str}（如不对请在上方修改）")

# pick files
sludge_name = pick_by_keywords(names, ["污泥"]) or pick_by_keywords(names, ["sludge"])
flow_name   = pick_by_keywords(names, ["流量"]) or pick_by_keywords(names, ["flow"])
turb_name   = pick_by_keywords(names, ["浊度"]) or pick_by_keywords(names, ["turb"])
ec_name     = pick_by_keywords(names, ["电导"]) or pick_by_keywords(names, ["EC"]) or pick_by_keywords(names, ["电导率"])
temp_name   = pick_by_keywords(names, ["温度"]) or pick_by_keywords(names, ["temp"])

missing=[]
for k,v in [("污泥",sludge_name),("流量",flow_name),("浊度",turb_name),("电导",ec_name),("温度",temp_name)]:
    if v is None: missing.append(k)

if missing:
    st.error("zip里没找到这些文件关键词：" + "、".join(missing) + "。请确认文件名包含对应中文关键字。")
    st.write("zip文件列表（前100个）：")
    st.write(names[:100])
    st.stop()

def read_inner(inner):
    with zipfile.ZipFile(uploaded, "r") as z:
        b=z.read(inner)
    return read_table_bytes(b, inner)

sludge=read_inner(sludge_name)
flow  =read_inner(flow_name)
turb  =read_inner(turb_name)
ec    =read_inner(ec_name)
temp  =read_inner(temp_name)

df = merge_site(date_str, sludge, flow, turb, ec, temp, tol_sec=int(tol_sec))

# quick NA diagnostics
na_rate = df.isna().mean().to_dict()
st.write("NA占比（越低越好）：", {k: f"{v*100:.1f}%" for k,v in na_rate.items()})

c0 = compute_c0(df)
ice_vol_val = ice_vol if ice_vol>0 else None
total_vol_val = total_vol if total_vol>0 else None

d, res = integrate(df, c0, ec_ice_thresh=ec_ice, ice_vol=ice_vol_val, total_vol=total_vol_val, expected_sec=expected_sec)

k = (report_mass / res["M_total_kg"]) if (report_mass>0 and res["M_total_kg"]>0) else None

st.subheader("结果摘要")
c1, c2, c3, c4 = st.columns(4)
c1.metric("C0(g/L)", f"{res['c0']:.4f}")
c2.metric("M_total(kg)", f"{res['M_total_kg']:.2f}")
c3.metric("M_ice(kg)", "—" if math.isnan(res["M_ice_kg"]) else f"{res['M_ice_kg']:.2f}")
c4.metric("k=报告/积分", "—" if k is None else f"{k:.2f}")

st.caption(f"数据质量：Q=0占比 {res['q0_ratio']*100:.1f}% ｜ 缺失占比(按{expected_sec}s) {res['missing_ratio']*100:.1f}%")

left, right = st.columns(2)
with left:
    st.plotly_chart(plot_timeseries(d), use_container_width=True)
with right:
    st.plotly_chart(plot_cumulative(d, res), use_container_width=True)

st.subheader("数据表（可下载）")
st.dataframe(d.head(200), use_container_width=True, height=360)

csv = d.to_csv(index=False).encode("utf-8-sig")
st.download_button("下载计算后数据CSV", data=csv, file_name="ice_flush_processed.csv", mime="text/csv")
